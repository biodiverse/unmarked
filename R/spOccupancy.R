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

  func <- quote(PGOcc)
  if(!is.null(umf@coordinates)){
    func <- quote(spPGOcc)
    data$coords <- umf@coordinates
  }

  forms <- split_formula(formula)
  
  #PGOcc(occ.formula = forms[[2]], det.formula = forms[[1]],
  #      data = data, ...)

  # Substitute version results in cleaner printed call
  cl <- substitute(
    spOccupancy::FUNC(occ.formula = OCCFORM, det.formula = DETFORM,
        data = data, ...),
    list(FUNC = func, OCCFORM = forms[[2]], DETFORM = forms[[1]]))
  
  # TODO: Allow outputing call and data instead of running?
  eval(cl)

})

setMethod("fit_spOcc", "unmarkedFrameOccuComm",
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
  
  M <- nrow(umf@ylist[[1]])
  J <- ncol(umf@ylist[[1]])
  S <- length(umf@ylist)
  yarr <- array(NA, c(S, M, J))
  rownames(yarr) <- names(umf@ylist)
  for (s in 1:S){
    yarr[s,,] <- umf@ylist[[s]]
  }

  data <- list(y = yarr, occ.covs = siteCovs(umf),
               det.covs = oc)

  func <- quote(msPGOcc)
  if(!is.null(umf@coordinates)){
    func <- quote(spMsPGOcc)
    data$coords <- umf@coordinates
  }

  vars <- all.vars(formula)
  if(any(vars %in% names(umf@speciesCovs))){
    stop("Species-level covariates not supported by spOccupancy", call.=FALSE)
  }

  forms <- split_formula(formula)
  
  #PGOcc(occ.formula = forms[[2]], det.formula = forms[[1]],
  #      data = data, ...)

  # Substitute version results in cleaner printed call
  cl <- substitute(
    spOccupancy::FUNC(occ.formula = OCCFORM, det.formula = DETFORM,
        data = data, ...),
    list(FUNC = func, OCCFORM = forms[[2]], DETFORM = forms[[1]]))
  
  # TODO: Allow outputing call and data instead of running?
  eval(cl)

})
