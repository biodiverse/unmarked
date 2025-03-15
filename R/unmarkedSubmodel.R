setClass("unmarkedSubmodel",
  slots = c(
    long_name = "character",
    short_name = "character",
    type = "character",
    data = "data.frame",
    formula = "formula",
    family = "character",
    link = "character",
    par_fixed = "numeric",
    vcov_fixed = "matrix",
    sigma = "data.frame",
    par_random = "data.frame",
    vcov_full = "matrix",
    auxiliary = "list"
  ),
  prototype = list(
    long_name = NA_character_,
    short_name = NA_character_,
    type = NA_character_,
    data = data.frame(),
    formula = ~1,
    family = NA_character_,
    link = NA_character_,
    par_fixed = numeric(0),
    vcov_fixed = matrix(0, 0, 0),
    sigma = data.frame(),
    par_random = data.frame(),
    vcov_full = matrix(0, 0, 0),
    auxiliary = list()
  )
)

setClass("unmarkedSubmodelState",
  contains = "unmarkedSubmodel")

setClass("unmarkedSubmodelDet",
  contains = "unmarkedSubmodel")

unmarkedSubmodelDet <- function(long_name, short_name, type, formula, data, 
                                family, link, auxiliary = list()){
  data <- clean_up_covs(data, drop_final = FALSE)$obs_covs
  data <- subset_covs(data, formula)
  new("unmarkedSubmodelDet",
    long_name = long_name, short_name = short_name, type = type,
    formula = formula, data = data, family = family, link = link,
    auxiliary = auxiliary)
}

unmarkedSubmodelState <- function(long_name, short_name, type, formula, data, 
                                  family, link, auxiliary = list()){
  data <- clean_up_covs(data, drop_final = FALSE)$site_covs
  data <- subset_covs(data, formula)
  new("unmarkedSubmodelState",
    long_name = long_name, short_name = short_name, type = type,
    formula = formula, data = data, family = family, link = link,
    auxiliary = auxiliary)
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


setMethod("model.matrix", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  form <- reformulas::nobars(object@formula)
  mf <- model.frame(form, object@data, na.action = na.pass)
  
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

#setGeneric("find_missing", function(object, ...){ 
#  standardGeneric("find_missing")
#})

#setMethod("find_missing", "unmarkedSubmodel",
#  function(object, ...){
#  mm <- model.matrix(object)
#  apply(mm, 1, function(x) any(is.na(x)))
#})

setGeneric("get_offset", function(object, ...){
  standardGeneric("get_offset")
})

setMethod("get_offset", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  mf <- model.frame(reformulas::nobars(object@formula), object@data, na.action=na.pass)
  if(!is.null(newdata)){
    X_terms <- stats::terms(mf)
    mf <- model.frame(X_terms, newdata, na.action=stats::na.pass)
  }
  off <- as.vector(stats::model.offset(mf))
  if(is.null(off)) off <- rep(0, nrow(mf))
  off
})

#setGeneric("get_parameter_idx", function(object, ...){
#  standardGeneric("get_parameter_idx")
#})

#setMethod("get_parameter_idx", "unmarkedSubmodelList",
#  function(object, ...){
#  idx <- lapply(object@submodels, function(x){
#    1:ncol(model.matrix(x))
#  })
#  for (i in 2:length(idx)){
#    idx[[i]] <- idx[[i]] + max(idx[[i-1]])
#  }
#  idx
#})

#setOldClass("sdreport")

setGeneric("add_fixed_estimates", function(object, tmb, sdr, ...){ 
  standardGeneric("add_fixed_estimates")
})

setMethod("add_fixed_estimates", "unmarkedSubmodel", function(object, tmb, sdr){   
  pars <- tmb$env$last.par.best
  par_fixed_idx <- names(pars) == paste0("beta_", object@type)
  par_fixed <- pars[par_fixed_idx]
  fixed_names <- colnames(model.matrix(object))
  names(par_fixed) <- fixed_names
  object@par_fixed <- par_fixed

  if(!is.null(sdr)){
    par_fixed_idx <- names(sdr$par.fixed) == paste0("beta_", object@type)
    vcov_fixed <- sdr$cov.fixed[par_fixed_idx, par_fixed_idx, drop = FALSE]
  } else {
    vcov_fixed <- matrix(NA, length(par_fixed), length(par_fixed))
  }
  colnames(vcov_fixed) <- rownames(vcov_fixed) <- fixed_names 
  object@vcov_fixed <- vcov_fixed

  object
})

setGeneric("has_random", function(object) standardGeneric("has_random"))

setMethod("has_random", "formula", function(object){
  length(reformulas::findbars(object)) > 0
})

setMethod("has_random", "unmarkedSubmodel", function(object){
  has_random(object@formula) 
})

setGeneric("add_random_estimates", function(object, tmb, sdr, ...){ 
  standardGeneric("add_random_estimates")
})

setMethod("add_random_estimates", "unmarkedSubmodel", function(object, tmb, sdr){

  if(!has_random(object)) return(object)

  # Random effects info
  re_info <- get_reTrms(object@formula, object@data)

  # Random effect variances (sigmas)
  pars <- tmb$env$last.par.best
  sig_fixed_idx <- names(pars) == paste0("lsigma_", object@type)
  log_sigma <- pars[sig_fixed_idx]

  if(!is.null(sdr)){
    sig_fixed_idx <- names(sdr$par.fixed) == paste0("lsigma_", object@type)
    log_sigma_var <- diag(sdr$cov.fixed)[sig_fixed_idx] 
  } else {
    log_sigma_var <- NA
  }

  rand_names <- random_effect_names(object)

  sigma_ests <- data.frame(Model = object@short_name,
                           rand_names$sigma,
                           log_sigma = log_sigma,
                           log_sigma_var = log_sigma_var)

  object@sigma <- sigma_ests 

  # Random effects
  ranef_idx <- names(pars) == paste0("b_", object@type)
  ranef_est <- pars[ranef_idx]

  ranef_est <- data.frame(Model = object@short_name,
                          rand_names$random_effects,
                          Estimate = ranef_est)

  object@par_random <- ranef_est

  # vcov matrix
  if(!is.null(sdr)){
    vcov_full <- solve(sdr$jointPrecision)
  } else {
    vcov_full <- matrix(NA, length(pars), length(pars)) 
  }  
  nms <- rownames(vcov_full)
  type_match <- grepl(paste0("_", object@type), nms, fixed=TRUE)
  not_sigma <- !grepl("lsigma_", nms, fixed = TRUE)
  idx <- type_match & not_sigma
  object@vcov_full <- vcov_full[idx, idx, drop = FALSE]

  object
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


setMethod("summary", "unmarkedSubmodel",
  function(object){

  ests <- object@par_fixed
  SEs <- SE(object)
  Z <- ests/SEs
  p <- 2*pnorm(abs(Z), lower.tail = FALSE)
  printRowNames <- !(length(ests) == 1 |
                       identical(names(ests), "(Intercept)") |
                       identical(names(ests), "1"))

  cat(object@long_name, " (", object@link, "-scale)", ":\n", sep="")

  if(has_random(object)){
    sigma_print <- object@sigma[,c("Groups", "Name")]
    sigma_print$Variance <- round(exp(object@sigma$log_sigma)^2, 3)
    sigma_print$`Std.Dev.` <- round(exp(object@sigma$log_sigma), 3)
    cat("Random effects:\n")
    print(sigma_print, row.names=FALSE)
    cat("\nFixed effects:\n")
  }

  outDF <- data.frame(Estimate = ests, SE = SEs, z = Z, "P(>|z|)" = p,
                      check.names = FALSE)
  print(outDF, row.names = printRowNames, digits = 3)
  cat("\n")
  invisible(outDF)
})

setMethod("show", "unmarkedSubmodel", function(object){
  summary(object)
})

#setMethod("AIC", "unmarkedFit", function(object){
#  2 * length(object@opt$par) - 2 * logLik(object)
#})

setGeneric("Z_matrix", function(object, ...){
  standardGeneric("Z_matrix")
})

setMethod("Z_matrix", "unmarkedSubmodel",
  function(object, newdata = NULL, ...){
  get_Z(object@formula, object@data, newdata)
})

setGeneric("get_TMB_pars", function(object, starts){
  standardGeneric("get_TMB_pars")
})

setMethod("get_TMB_pars", "unmarkedSubmodel", function(object, starts){

  out <- list()

  # Fixed pars
  fixed <- colnames(model.matrix(object))
  out$beta <- setNames(rep(0, length(fixed)), fixed)

  # Random effect names
  if(has_random(object)){
    re_info <- random_effect_names(object)

    # Sigma pars
    sig_names <- re_info$sigma
    sig_names <- apply(sig_names, 1, paste, collapse = "_")
    sig_names <- gsub("(Intercept)", "Int", sig_names, fixed = TRUE)
    out$lsigma <- setNames(rep(0, length(sig_names)), sig_names)

    # Random effect pars
    b_names <- apply(re_info$random_effects, 1, paste, collapse = "_")
    b_names <- gsub("(Intercept)", "Int", b_names, fixed = TRUE)
    out$b <- setNames(rep(0, length(b_names)), b_names)
  } else {
    out$lsigma <- numeric(0)
    out$b <- numeric(0)
  }
  
  names(out) <- paste0(names(out), "_", object@type)
  
  # Update start values
  if(!is.null(starts)){
    for (i in length(starts)){
      match_starts <- which(names(starts)[i] == names(out))
      if(length(match_starts) == 0) next
      if(length(starts[[i]]) != length(out[[match_starts]])){
        stop("starts element ", names(starts)[i], " should be length ",
             length(out[[match_starts]]), call.=FALSE)
      }
      out[[match_starts]] <- starts[[i]]
    }
  }

  out
})

setGeneric("get_TMB_random", function(object){
  standardGeneric("get_TMB_random")
})

setMethod("get_TMB_random", "unmarkedSubmodel", 
  function(object){
  if(!has_random(object)) return(NULL)
  paste0("b_", object@type)
})

setGeneric("get_TMB_data", function(object){
  standardGeneric("get_TMB_data")
})

setMethod("get_TMB_data", "list", function(object){
  object
})

setMethod("get_TMB_data", "unmarkedSubmodel", function(object){ 
  n_grouplevels <- get_nrandom(object@formula, object@data)
  out <- list(family = get_family_code(object@family),
              link = get_link_code(object@link),
              X = model.matrix(object), 
              offset = get_offset(object),
              Z = Z_matrix(object),
              n_group_vars = length(n_grouplevels[n_grouplevels > 0]),
              n_grouplevels = n_grouplevels)
  out <- c(out, object@auxiliary)
  names(out) <- paste0(names(out), "_", object@type)
  out
})

get_family_code <- function(family){
  switch(family, binomial = 0, poisson = 1)
}

get_link_code <- function(link){
  switch(link, logit = 0, exp = 1, cloglog = 2)
}

#find_gradient <- function(link){
#  switch(link,
#         logistic = "logistic.grad",
#         cloglog = "cloglog.grad",
#         exp = "exp")
#}

setMethod("predict", "unmarkedSubmodel",
  function(object, newdata = NULL, backTransform = TRUE, appendData = FALSE,
           level = 0.95, re.form = NULL, chunk_size = 70, ...){
  
  X <- model.matrix(object, newdata = newdata)
  coefs <- coef(object)
  off <- get_offset(object)
  cov_mat <- vcov(object)
  
  if(has_random(object) & is.null(re.form)){
    X <- cbind(X, Z_matrix(object, newdata = newdata))
    coefs <- coef(object, fixedOnly = FALSE)
    cov_mat <- vcov(object, fixedOnly = FALSE)
  }
  est <- as.vector(X %*% coefs)
  
  se <- lower <- upper <- rep(NA, length(est))
  if(!is.null(level)){
    ind <- rep(1:ceiling(nrow(X)/chunk_size), each=chunk_size, length.out=nrow(X)) 
    t_method <- base::t
    if(inherits(X, "Matrix")) t_method <- Matrix::t

    vars <- lapply(1:max(ind), function(i){
      x <- X[ind == i, , drop = FALSE]
      v <- as.matrix(x %*% cov_mat %*% t_method(x))
      diag(v)
    })
    vars <- do.call(c, vars)
    se <- sqrt(vars)
    z <- qnorm((1-level)/2, lower.tail = FALSE)
    lower <- est - z*se
    upper <- est + z*se
  }
  out <- data.frame(Predicted = est, SE = se, lower = lower, upper = upper,
                    row.names = NULL)

  if(backTransform){
    invlink <- get_invlink(object)
    out$Predicted <- invlink(out$Predicted)
    out$SE <- NA
    out$lower <- invlink(out$lower)
    out$upper <- invlink(out$upper)
  }

  if(appendData){
    out <- cbind(out, newdata)
  }

  out
})

get_invlink <- function(object){
  switch(object@link, 
         logit = plogis,
         log = exp,
         cloglog = function(x) 1-exp(-exp(x))
        )
}

setMethod("coef", "unmarkedSubmodel", 
  function(object, altNames = TRUE, fixedOnly = TRUE, ...){
 
  coefs <- object@par_fixed
  if(has_random(object) & !fixedOnly){
    rand <- object@par_random
    est <- object@par_random$Estimate
    nms <- paste(rand$Name, rand$Groups, rand$Levels, sep = "_")
    nms <- gsub("(Intercept)", "Int", nms, fixed=TRUE)
    names(est) <- nms
    coefs <- c(coefs, est)
  }

  if(altNames){
    names(coefs)[names(coefs) == "(Intercept)"] <- "Int"
    names(coefs) <- paste0(object@short_name, "(", names(coefs), ")")
  }

  coefs
})

setMethod("vcov", "unmarkedSubmodel", function(object, fixedOnly=TRUE, ...){
  v <- object@vcov_fixed
  if(!fixedOnly & nrow(object@vcov_full) > 1) v <- object@vcov_full
  rownames(v) <- colnames(v) <- names(coef(object, fixedOnly=fixedOnly))
  v
})

setMethod("SE", "unmarkedSubmodel", function(obj, fixedOnly = TRUE){
  v <- vcov(obj, fixedOnly = fixedOnly)
  sqrt(diag(v))
})

setMethod("confint", "unmarkedSubmodel", function(object, parm, level = 0.95,
                                                  fixedOnly = TRUE, ...){
  
  ests <- coef(object, fixedOnly = fixedOnly)
  ses <- SE(object, fixedOnly = fixedOnly)
  z <- qnorm((1-level)/2, lower.tail = FALSE)
  lower <- ests - z*ses
  upper <- ests + z*ses
  ci <- as.matrix(cbind(lower, upper))
  colnames(ci) <- c((1-level)/2, 1- (1-level)/2)

  if(missing(parm)) parm <- 1:length(ests)

  ci[parm,,drop=FALSE]
})

setMethod("sigma", "unmarkedSubmodel", function(object, level = 0.95, ...){
  if(!has_random(object)) stop("No random effects", call.=FALSE) 
  out <- object@sigma[,c("Model", "Groups", "Name")]
  out$Estimate <- exp(object@sigma$log_sigma)
  z <- qnorm((1-level)/2, lower.tail = FALSE)
  out$lower <- exp(out$Estimate - z * sqrt(object@sigma$log_sigma_var))
  out$upper <- exp(out$Estimate + z * sqrt(object@sigma$log_sigma_var))
  out
})


# submodel lists

setClass("unmarkedSubmodelList",
  slots = c(submodels = "list"),
  prototype = list(submodels = list())
)

unmarkedSubmodelList <- function(...){
  submodels <- list(...)
  names(submodels) <- sapply(submodels, function(x) x@type)
  new("unmarkedSubmodelList", submodels=submodels)
}

setGeneric("submodels", function(object) standardGeneric("submodels"))

setMethod("submodels", "unmarkedSubmodelList", function(object){
  object@submodels
})

setGeneric("submodels<-", function(object, value) standardGeneric("submodels<-"))

setMethod("submodels<-", "unmarkedSubmodelList", function(object, value){
  object@submodels <- value
  object
})

setMethod("get_TMB_random", "unmarkedSubmodelList",
  function(object){
  unlist(sapply(submodels(object), get_TMB_random))
})

setMethod("has_random", "unmarkedSubmodelList", function(object){
  sapply(submodels(object), has_random)
})

add_estimates <- function(object, tmb, se){
  sdr <- NULL
  if(se){
    sdr <- TMB::sdreport(tmb, getJointPrecision = TRUE)
  }
  submodels(object) <- lapply(submodels(object), add_fixed_estimates,
                             tmb = tmb, sdr = sdr)
  submodels(object) <- lapply(submodels(object), add_random_estimates,
                             tmb = tmb, sdr = sdr)
  object
}

setMethod("get_TMB_pars", "unmarkedSubmodelList", function(object, starts){
  unlist(lapply(unname(submodels(object)), get_TMB_pars, starts = starts), 
         recursive = FALSE)
})

setMethod("get_TMB_data", "unmarkedSubmodelList", function(object){
  out <- lapply(submodels(object), get_TMB_data)
  do.call(c, unname(out))
})

setGeneric("get_starts", function(object) standardGeneric("get_starts"))
setMethod("get_starts", "unmarkedSubmodelList", function(object){
  out <- lapply(submodels(object), get_TMB_pars, starts = NULL)
  rand <- sapply(submodels(object), get_TMB_random)
  out <- do.call(c, unname(out))
  out <- out[!names(out) %in% rand]
  out[sapply(out, length) > 0]
})



fit_TMB2 <- function(model, components, starts, se, method, ...){

  data <- get_TMB_data(components)
  params <- get_TMB_pars(submodelList(components), starts = starts)
  random <- get_TMB_random(submodelList(components))

  tmb_mod <- TMB::MakeADFun(data = c(model = model, data),
                            parameters = params,
                            random = random,
                            silent=TRUE,
                            DLL = "unmarked_TMBExports")

  opt <- optim(tmb_mod$par, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)

  submodelList(components) <- add_estimates(submodelList(components), tmb_mod, se)
  
  list(opt=opt, TMB=tmb_mod, components = components)
}

setAs("unmarkedSubmodel", "unmarkedEstimate", function(from){

  invlink <- switch(from@link, 
                    logit = "logistic", log = "exp", cloglog = "cloglog")

  grad <- switch(from@link, 
                 logit = "logistic.grad", log = "exp", cloglog = "cloglog.grad")

  new("unmarkedEstimate", name = from@long_name, short.name = from@short_name,
      estimates = coef(from, fixedOnly = FALSE),
      covMat = vcov(from, fixedOnly = FALSE),
      fixed = 1:length(coef(from)),
      invlink = invlink, invlinkGrad = grad, randomVarInfo = list())
})

setMethod("linearComb", c("unmarkedSubmodel", "matrixOrVector"),
  function(obj, coefficients, offset = NULL, re.form = NULL){
  est <- as(obj, "unmarkedEstimate")
  linearComb(est, coefficients = coefficients, offset = offset, re.form = re.form)
})

setMethod("backTransform", "unmarkedSubmodel", function(obj){
  est <- as(obj, "unmarkedEstimate")
  backTransform(est)
})
