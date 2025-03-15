setClass("unmarkedFit2",
  slots = c(
    components = "unmarkedModelComponents",
    opt = "list",
    TMB = "list",
    auxiliary = "list"
  ),
  prototype = list(
    opt = list(),
    TMB = list(),
    auxiliary = list()
  ),
  contains = "unmarkedFit"
)

setMethod("submodels", "unmarkedFit2", function(object){
  object@components@submodels@submodels
})

setMethod("summary_internal", "unmarkedFit2", function(object){
  cat("\nCall:\n")
  print(object@call)
  cat("\n")
  lapply(submodels(object), summary) 
  cat("AIC:", AIC(object),"\n\n")

  if(!identical(object@opt$convergence, 0L)){
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.", call.=FALSE)
  }

  # Check for potentially bad estimates
  if(is.null(object@call$se) || object@call$se){
    se <- SE(object)
    has_na <- any(is.na(se)) | any(is.nan(se))
    big_se <- any(abs(se) >= 5)
    if(has_na | big_se){
      warning("Large or missing SE values. Be very cautious using these results.", call.=FALSE)
    }
  }

})

setMethod("show", "unmarkedFit2", function(object){
  summary(object)
})

setMethod("[", "unmarkedFit2", function(x, i, j, drop) {
  submodels(x)[[i]]
})

setMethod("names", "unmarkedFit2", function(x){
  names(submodels(x))
})

setMethod("coef", "unmarkedFit2", 
  function(object, altNames = TRUE, fixedOnly = TRUE, ...){
  
  out <- lapply(submodels(object), coef,
                altNames = altNames, fixedOnly = fixedOnly)
  
  do.call(c, unname(out))
})

setMethod("vcov", "unmarkedFit2",
  function(object, fixedOnly=TRUE, method = c("hessian", "nonparboot"), ...){ 
  method <- match.arg(method)
  if(method == "nonparboot"){
    v <- covMatBS(object)
    nms <- names(coef(object, fixedOnly = fixedOnly))
    v <- v[rownames(v) %in% nms, colnames(v) %in% nms, drop = FALSE]
  } else {
    if(fixedOnly){
      sdr <- TMB::sdreport(object@TMB)
      is_sigma <- grepl("lsigma_", rownames(sdr$cov.fixed))
      v <- sdr$cov.fixed[!is_sigma, !is_sigma, drop = FALSE]
    } else {
      sdr <- TMB::sdreport(object@TMB, getJointPrecision = TRUE)
      is_sigma <- grepl("lsigma_", rownames(sdr$jointPrecision))
      v <- solve(sdr$jointPrecision)[!is_sigma, !is_sigma, drop = FALSE]
    }
    nms <- names(coef(object, fixedOnly = fixedOnly))
    rownames(v) <- colnames(v) <- nms
  }
  v
})

setMethod("SE", "unmarkedFit2", function(obj, fixedOnly = TRUE){
  out <- lapply(submodels(obj), SE, fixedOnly = fixedOnly)
  do.call(c, unname(out))
})

setMethod("predict", "unmarkedFit2",
  function(object, type, newdata = NULL, backTransform = TRUE, appendData = FALSE,
           level = 0.95, re.form = NULL, chunk_size = 70, ...){

  predict(object[type], newdata = newdata, backTransform = backTransform, 
    appendData = appendData, level = level, re.form = re.form,
    chunk_size = chunk_size, ...)
})

setMethod("confint", "unmarkedFit2",
  function(object, parm, level = 0.95, type, fixedOnly = TRUE, ...){

  if(!missing(type)){
    return(confint(object[type], parm = parm, level = level, fixedOnly = fixedOnly))
  }

  ci <- lapply(submodels(object), confint,
               parm = parm, level = level, fixedOnly = fixedOnly, ...)

  do.call(rbind, unname(ci))

})

setMethod("has_random", "unmarkedFit2", function(object){
  sapply(submodels(object), has_random)
})

setMethod("sigma", "unmarkedFit2", function(object, level = 0.95, ...){
  has_rand <- has_random(object)
  if(!any(has_rand)) stop("No random effects", call.=FALSE)
  sig <- lapply(submodels(object)[has_rand], sigma,
                level = level)
  do.call(rbind, unname(sig))
})

logLik.unmarkedFit2 <- function(object, ...){
  val <- -object@opt$value
  attr(val, "df") <- length(object@opt$par)
  class(val) <- "logLik"
  val
}

setMethod("linearComb", c("unmarkedFit2", "matrixOrVector"),
  function(obj, coefficients, type, offset = NULL, re.form = NULL){
  stopifnot(!missing(type))
  stopifnot(type %in% names(obj))
  estimate <- as(obj[type], "unmarkedEstimate")
  linearComb(estimate, coefficients, offset, re.form)
})

setMethod("backTransform", "unmarkedFit2", function(obj, type){
  est <- as(obj[type], "unmarkedEstimate")
  backTransform(est)
})
