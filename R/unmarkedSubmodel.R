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

setGeneric("has_random", function(object) standardGeneric("has_random"))

setMethod("has_random", "formula", function(object){
  length(reformulas::findbars(formula)) > 0
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
  mf <- model.frame(reformulas::nobars(object@formula), object@data, na.action=na.pass)
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
  idx
})

setMethod("engine_inputs", c("unmarkedSubmodelList", "missing"), function(object, object2){
  out <- lapply(submodels(object), engine_inputs)
  do.call(c, unname(out))
})

setMethod("engine_inputs_TMB", c("unmarkedSubmodelList", "missing"), function(object, object2){
  out <- lapply(submodels(object), engine_inputs_TMB)
  do.call(c, unname(out))
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


setGeneric("engine_inputs_R", function(object, object2) standardGeneric("engine_inputs_R"))

setMethod("engine_inputs_R", c("unmarkedResponse", "unmarkedSubmodelList"), 
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














fit_optim <- function(nll_fun, inputs, starts, method, se, ...){
  opt <- optim(starts, nll_fun, method = method, hessian = se, 
              inputs = inputs, ...)            
  cov_mat <- invertHessian(opt, length(opt$par), se)

  list(opt = opt, cov_mat = cov_mat,
       AIC = 2 * opt$value + 2 * length(starts)) #+ 2*nP*(nP + 1)/(M - nP - 1)
}








setGeneric("add_estimates", function(object, fit, ...) standardGeneric("add_estimates"))

# optim method (R/C engines)
setMethod("add_estimates", c("unmarkedSubmodelList", "list"),
  function(object, fit, ...){
  idx <- get_parameter_idx(object)
  for (i in names(submodels(object))){
    idx_rng <- idx[[i]][1]:idx[[i]][2] 
    object@estimates[[i]]@estimates <- fit$opt$par[idx_rng]
    v <- fit$cov_mat[idx_rng, idx_rng, drop = FALSE]
    object@estimates[[i]]@covMat <- v
  }
  object
})
