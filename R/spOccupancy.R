setGeneric("fit_spOcc", function(umf, ...) standardGeneric("fit_spOcc"))

# single-season single-species non-spatial occupancy
setMethod("fit_spOcc", "unmarkedFrameOccu",
          function(umf, formula, ...){
  
  if(!requireNamespace("spOccupancy")){
    stop("Install spOccupancy package", call.=FALSE)
  }
  oc <- list()
  if(!is.null(obsCovs(umf))){
    oc_add <- lapply(obsCovs(umf), function(x){
      matrix(x, nrow=numSites(umf), ncol=obsNum(umf), byrow=TRUE)
    })
    oc <- c(oc, oc_add)
  }
  if(!is.null(siteCovs(umf))){
    oc_add <- as.list(siteCovs(umf))
    oc <- c(oc, oc_add)
  }

  data <- list(y = getY(umf), occ.covs = siteCovs(umf),
               det.covs = oc)

  forms <- split_formula(formula)
  
  #PGOcc(occ.formula = forms[[2]], det.formula = forms[[1]],
  #      data = data, ...)

  # Substitute version results in cleaner printed call
  cl <- substitute(
    spOccupancy::PGOcc(occ.formula = OCCFORM, det.formula = DETFORM,
        data = data, ...),
    list(OCCFORM = forms[[2]], DETFORM = forms[[1]]))
  
  # TODO: Allow outputing call and data instead of running?
  eval(cl)

})
