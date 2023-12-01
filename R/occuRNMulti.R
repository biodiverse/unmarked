setClass("unmarkedFrameOccuRNMulti",
         contains="unmarkedFrameOccuMulti")

setClass("unmarkedFitOccuRNMulti",
         representation(
            detformulas = "list",
            stateformulas = "list"),
         contains = "unmarkedFit")

occuRNMulti <- function(detformulas, stateformulas, data, modelOccupancy,
                        K = rep(10, length(stateformulas)),
                        starts, method="BFGS", se = TRUE, threads=1, ...){
 
  stopifnot(all(names(data@ylist) == names(stateformulas)))

  S <- length(stateformulas)
  dep <- create_dep_matrix(stateformulas)

  if(missing(modelOccupancy)){
    modelOccupancy <- rep(0, S)
    names(modelOccupancy) <- names(stateformulas)
  } else {
    stopifnot(identical(names(modelOccupancy), names(stateformulas)))
    if(is.logical(modelOccupancy)){
      modelOccupancy <- modelOccupancy*1
    }
    not_allowed_occ <- colSums(dep)
    sp_mismatch <- modelOccupancy & not_allowed_occ
    if(any(sp_mismatch)){
      stop(paste("Species", 
                 paste(names(modelOccupancy)[sp_mismatch], collapse=", "),
                 "cannot use occupancy models."), call.=FALSE)
    }
  }
  
  data <- as(data, "unmarkedFrameOccuRNMulti")
  gd <- getDesign(data, detformulas, stateformulas)

  Kmin <- sapply(gd$ylist, function(x) apply(x, 1, max))

  nP <- max(gd$det_ind)
  
  has_occ <- sum(modelOccupancy) > 0
  if(has_occ){
    sp_occ <- names(modelOccupancy)[modelOccupancy == 1]
    sp_abun <- names(modelOccupancy)[modelOccupancy == 0]
    occ_ind <- gd$state_ind[grepl(paste(paste0("^",sp_occ), collapse="|"),
                                  rownames(gd$state_ind)),,drop=FALSE]
    occ_rng <- min(occ_ind):max(occ_ind)
    abun_ind <- gd$state_ind[grepl(paste(paste0("^",sp_abun), collapse="|"),
                                  rownames(gd$state_ind)),,drop=FALSE]
    abun_rng <- min(abun_ind):max(abun_ind)
  } else {
    state_rng <- min(gd$state_ind):max(gd$state_ind)
  }

  det_rng <- min(gd$det_ind):max(gd$det_ind)
  
  if(missing(starts)) starts <- rep(0, nP)
  names(starts) <- gd$par_names

  nll_C <- function(pars){
    nll_occuRNMulti(pars, gd$state_ind-1, gd$det_ind-1, S, modelOccupancy,
                    gd$ylist, gd$state_dm, gd$det_dm, dep, K, Kmin, threads)
  }

  fm <- optim(starts, nll_C, method=method, hessian=se, ...)
  covMat <- invertHessian(fm, nP, se)

  fmAIC <- 2 * fm$value + 2 * nP
  #----------------------------------------------------------------------------

  #Format output---------------------------------------------------------------
  ests <- fm$par
  names(ests) <- gd$par_names

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates = ests[det_rng],
                          covMat = as.matrix(covMat[det_rng, det_rng]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

  if(has_occ){
    lam <- unmarkedEstimate(name = "Abundance", short.name = "lam",
                            estimates = ests[abun_rng],
                            covMat = as.matrix(covMat[abun_rng, abun_rng]),
                            invlink = "exp",
                            invlinkGrad = "exp")
    psi <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                          estimates = ests[occ_rng],
                          covMat = as.matrix(covMat[occ_rng, occ_rng]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

    estimateList <- unmarkedEstimateList(list(lam=lam, psi=psi, det=det))
  } else {
    lam <- unmarkedEstimate(name = "Abundance", short.name = "lam",
                            estimates = ests[state_rng],
                            covMat = as.matrix(covMat[state_rng, state_rng]),
                            invlink = "exp",
                            invlinkGrad = "exp")
    estimateList <- unmarkedEstimateList(list(lam=lam, det=det))
  }

  umfit <- new("unmarkedFitOccuRNMulti", fitType = "occuRNMulti", call = match.call(),
                detformulas = detformulas, stateformulas = stateformulas,
                formula = ~1, data = data,
                #sitesRemoved = designMats$removed.sites,
                estimates = estimateList, AIC = fmAIC, opt = fm,
                negLogLike = fm$value, nllFun = nll_C)

  return(umfit)
}

setMethod("getDesign", "unmarkedFrameOccuRNMulti",
  function(umf, detformulas, stateformulas){
  M <- nrow(umf@ylist[[1]])
  J <- ncol(umf@ylist[[1]])
  
  ylist <- umf@ylist

  # Site covs
  sc <- umf@siteCovs
  if(is.null(sc)) sc <- data.frame(dummy_=rep(1, M))

  mm <- lapply(stateformulas, function(s){
    if(!is.list(s)){
      return(model.matrix(s, sc))
    }
    lapply(s, function(x) model.matrix(x, sc))
  })

  # Obs covs
  sc_long_ind <- rep(1:M, each=J)
  sc_long <- sc[sc_long_ind,]
  oc <- cbind(sc_long, umf@obsCovs)

  vv <- lapply(detformulas, function(x){
    model.matrix(x, oc)
  })

  # Indices
  state_n <- unlist(lapply(mm, function(x){
    if(is.list(x)){
      lapply(x, function(y) ncol(y))
    } else {
      ncol(x)
    }
  }))
  state_end <- cumsum(state_n)
  state_start <- state_end - state_n + 1
  state_ind <- cbind(state_start, state_end)

  state_prefix <- rep(rownames(state_ind), state_n)
  state_prefix <- gsub(".", ":", state_prefix, fixed=TRUE)
  state_prefix <- paste0("[",state_prefix,"]")
  
  state_par <- unlist(lapply(mm, function(x){
    if(is.matrix(x)) return(colnames(x))
    lapply(x, function(y){
      return(colnames(y))
    })
  }))
  stopifnot(length(state_prefix) == length(state_par))
  state_par <- paste(state_prefix, state_par)
  
  det_n <- sapply(vv, ncol)
  det_end <- cumsum(det_n) + max(state_end)
  det_start <- det_end - det_n + 1
  det_ind <- cbind(det_end, det_start)

  det_rng <- min(det_ind):max(det_ind)

  det_prefix <- paste0("[",rep(rownames(det_ind), det_n),"]")
  det_par <- unlist(lapply(vv, colnames))
  stopifnot(length(det_prefix) == length(det_par))
  det_par <- paste(det_prefix, det_par)
  par_names <- c(state_par, det_par)

  list(par_names = par_names, state_ind = state_ind, det_ind = det_ind,
              ylist = ylist, state_dm = mm, det_dm = vv)
})

create_dep_matrix <- function(stateformulas){
  S <- length(stateformulas)
  dep <- matrix(0, S, S)
  rownames(dep) <- colnames(dep) <- names(stateformulas)
  allowed_species <- character(0)

  # Species 1
  stopifnot(is.call(stateformulas[[1]]))
  allowed_species <- c(allowed_species, names(stateformulas)[1])

# Species 2
  if(is.list(stateformulas[[2]])){
    dep_sp <- names(stateformulas[[2]])[2:length(stateformulas[[2]])]
    stopifnot(all(dep_sp %in% allowed_species))
    for (i in dep_sp){
      dep[names(stateformulas)[2], dep_sp] <- 1
    }
  }
  allowed_species <- c(allowed_species, names(stateformulas)[2])

  # Species 3
  if(S > 2){
    if(is.list(stateformulas[[3]])){
      dep_sp <- names(stateformulas[[3]])[2:length(stateformulas[[3]])]
      stopifnot(all(dep_sp %in% allowed_species))
      for (i in dep_sp){
        dep[names(stateformulas)[3], dep_sp] <- 1
      }
    }
    allowed_species <- c(allowed_species, names(stateformulas)[3])
  }

  dep
}


