setClass("unmarkedFrameOccuComm",
         representation(ylist = "list", speciesCovs="optionalList"),
         contains = "unmarkedFrame")

setClass("unmarkedFitOccuComm", contains="unmarkedFitOccu")


unmarkedFrameOccuComm <- function(y, siteCovs=NULL, obsCovs=NULL, speciesCovs=NULL){  
  # Check species covs dimensions
  if(!is.null(speciesCovs)){
    stopifnot(is.list(speciesCovs))
    M <- nrow(y[[1]])
    J <- ncol(y[[1]])
    S <- length(y)
    correct_dims <- sapply(speciesCovs, function(x){
      identical(dim(x), c(M, S)) | identical(dim(x), c(M, J, S)) 
    })
    bad_dims <- names(speciesCovs)[!correct_dims]
    if(length(bad_dims > 0)){
      stop(paste0("Species covariate(s) ", paste(bad_dims, collapse=", "), " have incorrect dimensions"),
           call.=FALSE)
    }
  }

  obsCovs <- covsToDF(obsCovs, "obsCovs", ncol(y[[1]]), nrow(y[[1]]))
  new("unmarkedFrameOccuComm", y=y[[1]], ylist = y, siteCovs=siteCovs, 
      obsCovs=obsCovs, speciesCovs=speciesCovs)
}

process_multispecies_umf <- function(umf, interact_covs){
  ylist <- umf@ylist
  M <- nrow(ylist[[1]])
  J <- ncol(ylist[[1]])
  S <- length(ylist)
  sp <- factor(names(ylist), levels=names(ylist))

  y <- do.call(rbind, ylist)

  sc <- siteCovs(umf)
  sp_rep <- rep(sp, each = M)
  if(is.null(sc)){
    sc <- data.frame(species=sp_rep)
  } else {
    sc <- sc[rep(1:M, S),,drop=FALSE]
    sc[["species"]] <- sp_rep

    # Handle random effects
    interact_covs <- names(sc)[names(sc) %in% interact_covs]
    for (i in interact_covs){
      stopifnot(is.factor(sc[[i]]) | is.character(sc[[i]]))
      sc[[paste0(i, "_species")]] <- factor(paste0(sc[[i]], "_", sc[["species"]]))  
    }
  }

  oc <- obsCovs(umf)
  #sp_rep <- rep(sp, each=M*J)
  if(!is.null(oc)){
    oc <- oc[rep(1:(M*J), S),,drop=FALSE]

    # Handle random effects
    interact_covs <- names(oc)[names(oc) %in% interact_covs]
    for (i in interact_covs){
      stopifnot(is.factor(oc[[i]]) | is.character(oc[[i]]))
      sp_rep <- rep(sc$species, each=M*J)
      oc[[paste0(i, "_species")]] <- factor(paste0(oc[[i]], "_", oc[["species"]]))    
    }
  }

  # Add species level covs
  spc <- umf@speciesCovs
  if(!is.null(spc)){
    # Site level species covs
    spc_site <- sapply(spc, function(x) identical(dim(x), c(M, S))) 
    spc_site <- spc[spc_site]
    if(length(spc_site) > 0){
      spc_site <- as.data.frame(lapply(spc_site, as.vector))
      if(is.null(sc)){
        sc <- spc_site
      } else {
        sc <- cbind(sc, spc_site)
      }
    }

    spc_obs <- sapply(spc, function(x) identical(dim(x), c(M, J, S)))
    spc_obs <- spc[spc_obs]
    if(length(spc_obs) > 0){
      spc_obs <- lapply(spc_obs, function(x){
        x <- aperm(x, c(2, 1, 3))
        as.vector(x)
      })
      spc_obs <- as.data.frame(spc_obs)
      if(is.null(oc)){
        oc <- spc_obs
      } else {
        oc <- cbind(oc, spc_obs)
      }
    }
  }

  unmarkedFrameOccu(y = y, siteCovs=sc, obsCovs=oc)
}

multispeciesFormula <- function(form){
  sf <- split_formula(form)

  covs_add_species <- character(0)

  # Det formula
  det_nobar <- reformulas::nobars(sf[[1]])
  bars <- ifelse(det_nobar == ~1, "|", "||")
  rand <- paste0("+ (", safeDeparse(det_nobar[[2]]), " ", bars, " species)")
  
  barexp <- reformulas::findbars(sf[[1]])
  if(!is.null(barexp)){
    lhs <- lapply(barexp, function(x) x[[2]])
    check_lhs <- sapply(lhs, function(x) x == 1)
    if(any(!check_lhs)) stop("Only random intercepts allowed", call.=FALSE)
    new_rhs <- lapply(barexp, function(x){
      rhs <- safeDeparse(x[[3]])
      covs_add_species <<- c(covs_add_species, rhs)
      new_rhs <- paste0(rhs, "_species")
      x[[3]] <- str2lang(new_rhs)
      x
    })
    rand2 <- lapply(new_rhs, function(x){
      paste0("(", safeDeparse(x), ")")
    })
    rand <- paste(rand, "+", paste(rand2, collapse="+"))
  }

  sf[[1]] <- paste(safeDeparse(det_nobar), rand)

  # State formula
  state_nobar <- reformulas::nobars(sf[[2]])
  bars <- ifelse(state_nobar == ~1, "|", "||")
  rand <- paste0("+ (", safeDeparse(state_nobar[[2]]), " ", bars, " species)")

  barexp <- reformulas::findbars(sf[[2]])
  if(!is.null(barexp)){
    lhs <- lapply(barexp, function(x) x[[2]])
    check_lhs <- sapply(lhs, function(x) x == 1)
    if(any(!check_lhs)) stop("Only random intercepts allowed", call.=FALSE)
    new_rhs <- lapply(barexp, function(x){
      rhs <- safeDeparse(x[[3]])
      covs_add_species <<- c(covs_add_species, rhs)
      new_rhs <- paste0(rhs, "_species")
      x[[3]] <- str2lang(new_rhs)
      x
    })
    rand2 <- lapply(new_rhs, function(x){
      paste0("(", safeDeparse(x), ")")
    })
    rand <- paste(rand, "+", paste(rand2, collapse="+"))
  }

  sf[[2]] <- paste(safeDeparse(state_nobar), rand)
  
  list(formula = stats::as.formula(str2lang(paste0(sf[[1]], sf[[2]]))),
       covs = unique(covs_add_species))
}

safeDeparse <- function(inp) {
  out <- deparse(inp)
  paste(sapply(out, trimws), collapse=" ")
}

occuComm <- function(formula, data, ...){
  newform <- multispeciesFormula(formula)
  newumf <- process_multispecies_umf(data, newform$covs)
  out <- occu(newform$formula, data=newumf, ...)
  out@call <- match.call()
  out@formula <- formula
  out@data <- data
  out <- as(out, "unmarkedFitOccuComm")
  out
}

setMethod("ranef_internal", "unmarkedFitOccuComm", function(object){
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")
  r <- ranef(new_object)
  inds <- split(1:nrow(r@post), rep(1:S, each=M))
  names(inds) <- names(object@data@ylist)
  lapply(inds, function(x){ 
         out <- r
         out@post <- out@post[x,,,drop=FALSE]
         out
  })
})

setGeneric("richness", function(object, ...) standardGeneric("richness"))
setMethod("richness", "unmarkedFitOccuComm", 
          function(object, nsims=100, posterior=FALSE){
  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")
  r <- ranef(new_object)
  post <- posteriorSamples(r, nsims=nsims)@samples
  post <- array(post, c(M, S, 1, nsims))
  rich <- apply(post, c(1,3,4), sum, na.rm=TRUE)
  if(!posterior){
    rich <- drop(rich)
    return(apply(rich, 1, mean, na.rm=TRUE))
  }
  new("unmarkedPostSamples", numSites=M, numPrimary=1, nsims=nsims,
      samples=rich)
})

setMethod("predict_internal", "unmarkedFitOccuComm",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){

  S <- length(object@data@ylist)
  M <- numSites(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")

  pr <- predict(new_object, type=type, newdata=newdata, backTransform=backTransform,
                na.rm=na.rm, appendData=appendData, level=level, re.form=re.form, ...)

  inds <- split(1:nrow(pr), rep(1:length(object@data@ylist), each=numSites(object@data)))
  names(inds) <- names(object@data@ylist)
  lapply(inds, function(x){ 
         out <- pr[x,,drop=FALSE]
         rownames(out) <- NULL
         out
  })
})
