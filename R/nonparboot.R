setGeneric("nonparboot", function(object, ...) standardGeneric("nonparboot"))

# Overall nonparboot method; this is the one exported
setMethod("nonparboot", "unmarkedFit", function(object, B = 1, keepOldSamples = TRUE, ...){
  bsType <- list(...)$bsType
  if(!is.null(bsType) && bsType == "both"){
    warning("Only bsType = 'site' is supported")
  }
  stopifnot(B > 0)
  # Calls fit-specific internal method
  nonparboot_internal(object, B=B, keepOldSamples=keepOldSamples)
})

# Fit-specific internal methods
setGeneric("nonparboot_internal", function(object, B, keepOldSamples){
  standardGeneric("nonparboot_internal")
})

# Select/insert bootstrap samples
setGeneric("bootstrapSamples", function(object) standardGeneric("bootstrapSamples"))

setMethod("bootstrapSamples", "unmarkedFit", function(object){
  object@bootstrapSamples
})

setMethod("bootstrapSamples", "unmarkedFit2", function(object){
  object@auxiliary$bootstrapSamples
})

setGeneric("bootstrapSamples<-", function(object, value){
  standardGeneric("bootstrapSamples<-")
})

setMethod("bootstrapSamples<-", "unmarkedFit", function(object, value){
  object@bootstrapSamples <- value
  object
})

setMethod("bootstrapSamples<-", "unmarkedFit2", function(object, value){
  object@auxiliary$bootstrapSamples <- value
  object
})

# Select/insert bootstrap varcov matrix
setGeneric("covMatBS", function(object) standardGeneric("covMatBS"))

setMethod("covMatBS", "unmarkedFit", function(object){
  object@covMatBS
})

setMethod("covMatBS", "unmarkedFit2", function(object){
  object@auxiliary$covMatBS
})

setGeneric("covMatBS<-", function(object, value){
  standardGeneric("covMatBS<-")
})

setMethod("covMatBS<-", "unmarkedFit", function(object, value){
  object@covMatBS <- value

  short_names <- regmatches(colnames(value), regexpr("^[^\\(]+", colnames(value)))
  
  for (i in 1:length(object@estimates@estimates)){
    match_name <- object@estimates@estimates[[i]]@short.name
    inds <- short_names == match_name
    new_v <- value[inds, inds, drop=FALSE]
    object@estimates@estimates[[i]]@covMatBS <- new_v
  }

  object
})

setMethod("covMatBS<-", "unmarkedFit2", function(object, value){
  object@auxiliary$covMatBS <- value
  object
})

setMethod("nonparboot_internal", "unmarkedFit", 
          function(object, B, keepOldSamples){

  M <- numSites(object@data)
  
  # Note: no sites are removed. I think this is more accurate as it
  # results in comparable sample sizes to original data.
  boot_iter <- lapply2(1:B, function(i){
    finish <- FALSE
    while(!finish){
      sites <- sort(sample(1:M, M, replace=TRUE))
      new_data <- getData(object)[sites,]
      ran <- TRUE
      tryCatch(fit <- nonparboot_update(object, data = new_data),
               error = function(e) ran <<- FALSE)
      conv <- object@opt$convergence == 0
      if(ran & conv) finish <- TRUE
    }
    fit
  })

  if(!keepOldSamples) bootstrapSamples(object) <- NULL
  bootstrapSamples(object) <- c(bootstrapSamples(object), boot_iter)
   
  coefs <- t(sapply(bootstrapSamples(object), coef, fixedOnly = FALSE))
  v <- stats::cov(coefs)
  covMatBS(object) <- v

  object
})

# Fit-specific update methods (to handle occuPEN)
setGeneric("nonparboot_update", function(object, data){
  standardGeneric("nonparboot_update")
})

setMethod("nonparboot_update", "unmarkedFit", function(object, data){
  update(object, data = data, se = FALSE)
})

setMethod("nonparboot_update", "unmarkedFitOccuPEN", function(object, data){
  update(object, data = data)
})

setMethod("nonparboot_update", "unmarkedFitOccuPEN_CV", function(object, data){
  if(object@pen.type == "MPLE"){
	  MPLElambda <- computeMPLElambda(object@formula, data)
    out <- update(object, data = data, lambda = MPLElambda)
  } else {
    out <- update(object, data = data)
  }
  out
})

# Colext has some additional project stuff
setMethod("nonparboot_internal", "unmarkedFitColExt", 
          function(object, B, keepOldSamples){

  # Call base method
  object <- methods::callNextMethod(object, B=B, keepOldSamples=keepOldSamples)
  
  # Add smoothed/projected stuff
  smoothed.occ <- t(sapply(object@bootstrapSamples,
                    function(x) x@smoothed.mean[1,]))
  smoothed.unocc <- t(sapply(object@bootstrapSamples,
                      function(x) x@smoothed.mean[2,]))
  object@smoothed.mean.bsse <- rbind(sqrt(diag(cov(smoothed.occ))),
                                     sqrt(diag(cov(smoothed.unocc))))
  projected.occ <- t(sapply(object@bootstrapSamples,
                     function(x) x@projected.mean[1,]))
  projected.unocc <- t(sapply(object@bootstrapSamples,
                    function(x) x@projected.mean[2,]))
  object@projected.mean.bsse <- rbind(sqrt(diag(cov(projected.occ))),
                                      sqrt(diag(cov(projected.unocc))))
  object
})
