setClass("unmarkedFrameOccuRNMulti",
         contains="unmarkedFrameOccuMulti")

setClass("unmarkedFitOccuRNMulti",
         representation(
            detformulas = "list",
            stateformulas = "list",
            modelOccupancy = "numeric"),
         contains = "unmarkedFit")

occuRNMulti <- function(detformulas, stateformulas, data, modelOccupancy,
                        K = rep(10, length(stateformulas)),
                        starts, method="BFGS", se = TRUE, threads=1, ...){
 
  stopifnot(all(names(data@ylist) == names(stateformulas)))

  S <- length(stateformulas)
  dep <- create_dep_matrix(stateformulas)

  modocc <- rep(0, S)
  names(modocc) <- names(data@ylist)
  if(!missing(modelOccupancy)){
    stopifnot(all(modelOccupancy %in% names(modocc)))
    modocc[modelOccupancy] <- 1
    not_allowed_occ <- colSums(dep)
    sp_mismatch <- modocc & not_allowed_occ
    if(any(sp_mismatch)){
      stop(paste("Species", 
                 paste(names(modocc)[sp_mismatch], collapse=", "),
                 "cannot use occupancy model."), call.=FALSE)
    }
  }

  data <- as(data, "unmarkedFrameOccuRNMulti")
  gd <- getDesign(data, detformulas, stateformulas)

  Kmin <- sapply(gd$ylist, function(x) apply(x, 1, max))

  nP <- max(gd$det_ind)
  
  has_occ <- sum(modocc) > 0
  if(has_occ){
    sp_occ <- names(modocc)[modocc == 1]
    sp_abun <- names(modocc)[modocc == 0]
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
    nll_occuRNMulti(pars, gd$state_ind-1, gd$det_ind-1, S, modocc,
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
                modelOccupancy = modocc,
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
  if(is.null(umf@obsCovs)){
    oc <- sc_long
  } else {
    oc <- cbind(sc_long, umf@obsCovs)
  }

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
  det_ind <- cbind(det_start, det_end)

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

setMethod("predict", "unmarkedFitOccuRNMulti",
          function(object, type, species, newdata,
                   level = 0.95, nsims = 100, ...){

  stopifnot(type %in% c("state", "det"))
  if(!missing(species)){
    stopifnot(species %in% names(object@data@ylist))
  }
  if(type == "state"){
    samps <- MASS::mvrnorm(nsims, mu=coef(object), Sigma = vcov(object))

    depmat <- unmarked:::create_dep_matrix(object@stateformulas)

    top <- apply(depmat, 1, function(x) sum(x) == 0)
    sp_top <- colnames(depmat)[top]
    lev2 <- apply(depmat[!top,!top,drop=FALSE], 1, function(x) sum(x) == 0)
    sp_2nd <- names(lev2)[lev2]
    sp_3rd <- names(lev2)[!lev2]

    if(missing(newdata) || is.null(newdata)){
      
      X <- lapply(object@stateformulas, function(x){
        if(is.list(x)){
          lapply(x, function(z) model.matrix(z, object@data@siteCovs))
        } else {
          model.matrix(x, object@data@siteCovs)
        }
      })
      nr <- nrow(object@data@siteCovs)
    } else {
      X <- lapply(object@stateformulas, function(x){
        if(is.list(x)){
          lapply(x, function(z) unmarked:::make_mod_matrix(z, object@data@siteCovs, newdata)$X)
        } else {
          unmarked:::make_mod_matrix(x, object@data@siteCovs, newdata)$X
        }
      })
      nr <- nrow(newdata)
    }
    
    coef_est <- matrix(coef(object), nrow=1)
    colnames(coef_est) <- names(coef(object))
    point_est <- calc_dependent_response(coef_est, X, object, nr, sp_top, sp_2nd, sp_3rd, ilink)

    message('Bootstrapping confidence intervals with ',nsims,' samples')
    post <- calc_dependent_response(samps, X, object, nr, sp_top, sp_2nd, sp_3rd, ilink)

    # Summarize bootstrap
    cis <- lapply(post, function(x){
      data.frame(
        SE = apply(x, 1, sd, na.rm=TRUE),
        upper = apply(x, 1, quantile, 0.025, na.rm=TRUE),
        lower = apply(x, 1, quantile, 0.975, na.rm=TRUE)
      )
    })

    out <- mapply(function(x, y){
      cbind(Predicted = x, y)
    }, point_est, cis, SIMPLIFY=FALSE)

    if(!missing(species)){
      out <- out[[species]]
    }
    return(out)
  } else if(type == "det"){
    stop("det not supported yet")
  }
})

# Function to extract correct parts of state coefs vector
get_species_b <- function(b, object, sp){
  b <- b[!grepl("^p\\(", names(b))] # remove detection coefficients
  nm <- names(b)
  allsp <- names(object@data@ylist)
  ind <- which(allsp == sp)
  out <- list()
  other <- c(1:length(allsp))[-ind]
  match1 <- grepl(paste0("[", allsp[ind], "]"), nm, fixed=TRUE)
  out[[1]] <- b[match1]

  match2 <- grepl(paste0("[", sp, ":", allsp[other[1]],"]"), nm, fixed=TRUE)
  if(sum(match2 > 0)){
    out[[allsp[other[1]]]] <- b[match2]
  }
  match3 <- grepl(paste0("[", sp, ":", allsp[other[2]],"]"), nm, fixed=TRUE)
  if(sum(match3 > 0)){
    out[[allsp[other[2]]]] <- b[match3]
  }
  out
}

calc_dependent_response <- function(samps, X, object, nr, sp_top, sp_2nd, sp_3rd, ilink){
  nsims <- nrow(samps)

  all_species <- names(object@data@ylist)

  post <- lapply(1:length(all_species), function(x) matrix(NA, nr, nsims))
  names(post) <- all_species

  for (k in 1:nsims){

    b <- samps[k,]
    ests <- list()

    # Top-level species that don't depend on other species abundance
    for (i in sp_top){
      beta <- get_species_b(b, object, i)[[1]]
      ests[[i]] <- exp(X[[i]] %*% beta)
      post[[i]][,k] <- ests[[i]]
    }
    # 2nd-level species that depend only on top-level species
    for (i in sp_2nd){
      if(object@modelOccupancy[i]){
        ilink <- plogis
      } else {
        ilink <- eval(str2lang(object@estimates@estimates$lam@invlink))
      }
      Xsub <- X[[i]]
      beta <- get_species_b(b, object, i)
      lp <- Xsub[[1]] %*% beta[[1]]
      for (j in 2:length(Xsub)){
        other_sp <- names(Xsub)[j]
        add <- Xsub[[other_sp]] %*% beta[[other_sp]] * ests[[other_sp]]
        lp <- lp + add
      }
      ests[[i]] <- ilink(lp)
      post[[i]][,k] <- ests[[i]]
    }
    # 3rd-level species that depend on at least one 2nd-level species
    for (i in sp_3rd){
      if(object@modelOccupancy[i]){
        ilink <- plogis
      } else {
        ilink <- eval(str2lang(object@estimates@estimates$lam@invlink))
      }
      Xsub <- X[[i]]
      beta <- get_species_b(b, object, i)
      lp <- Xsub[[1]] %*% beta[[1]]
      for (j in 2:length(Xsub)){
        other_sp <- names(Xsub)[j]
        add <- Xsub[[other_sp]] %*% beta[[other_sp]] * ests[[other_sp]]
        lp <- lp + add
      }
      ests[[i]] <- ilink(lp)
      post[[i]][,k] <- ests[[i]]
    }
  }
  post
}
