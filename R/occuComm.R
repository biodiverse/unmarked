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

  if(is.null(names(y))){
    names(y) <- paste0("sp", sprintf("%02d", 1:S))
  }

  obsCovs <- covsToDF(obsCovs, "obsCovs", ncol(y[[1]]), nrow(y[[1]]))
  new("unmarkedFrameOccuComm", y=y[[1]], ylist = y, siteCovs=siteCovs, 
      obsCovs=obsCovs, speciesCovs=speciesCovs, obsToY = diag(J))
}

# Unmarked Frame methods
setMethod("show", "unmarkedFrameOccuComm", function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    cat("Only showing observation matrix for species 1.\n")
    print(df)
})

setMethod("summary", "unmarkedFrameOccuComm", function(object,...) {
    cat("unmarkedFrameOccuComm Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat(length(object@ylist), "species\n")
    cat("Maximum number of observations per site:",obsNum(object),"\n")
    mean.obs <- mean(rowSums(!is.na(object@ylist[[1]])))
    cat("Mean number of observations per site:",round(mean.obs,2),"\n")
    s.one <- sapply(object@ylist,function(x) sum(rowSums(x,na.rm=TRUE)>0))
    cat("Sites with at least one detection, quantiles by species:\n")
    print(quantile(s.one))
    if(!is.null(object@siteCovs)) {
        cat("\nSite-level covariates:\n")
        print(summary(object@siteCovs))
    }
    if(!is.null(object@obsCovs)) {
        cat("\nObservation-level covariates:\n")
        print(summary(object@obsCovs))
    }
    if(!is.null(object@speciesCovs)){
        cat("\nSpecies-level covariates:\n")
        str(object@speciesCovs)
    }
})

setMethod("plot", c(x="unmarkedFrameOccuComm", y="missing"),
	function (x, y, colorkey, ylab="Site", xlab="Observation", ...)
{
    y <- getY(x)
    ym <- max(y, na.rm=TRUE)
    M <- nrow(y)
    J <- ncol(y)
    S <- length(x@ylist)
    y <- as.data.frame(do.call(rbind,x@ylist))
    colnames(y) <- paste("obs",1:J)
    y$site <- rep(1:M,S)
    y$species <- as.factor(rep(names(x@ylist),each=M))
    y2 <- reshape(y, idvar=c("site", "species"), varying=list(1:obsNum(x)),
              v.names="value", direction="long")
    y2$variable <- factor(paste("obs", y2$time), levels=colnames(y))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym),
            at=(0:ym)+0.5))
    levelplot(value ~ variable*site | species, y2,
        scales=list(relation="free", x=list(labels=1:J)),
        colorkey=colorkey, strip=T, xlab=xlab, ylab=ylab,
        labels=names(x@ylist), ...)
})

#[ Methods for community occupancy frames
setMethod("[", c("unmarkedFrameOccuComm", "numeric", "missing", "missing"),
    function(x, i){
  if(length(i) == 0) return(x)
  M <- numSites(x)
  J <- obsNum(x)
  S <- length(x@ylist)

  ylist <- lapply(x@ylist,function(x) x[i,,drop=F])
 
  siteCovs <- siteCovs(x)
  if (!is.null(siteCovs)) {
    siteCovs <- siteCovs(x)[i, , drop = FALSE]
  }

  obsCovs <- obsCovs(x)
  if (!is.null(obsCovs)) {
    .site <- rep(1:M, each = J)
    oc <- lapply(i, function(ind){
      obsCovs[.site==ind,,drop=FALSE]
    })
    obsCovs <- do.call(rbind, oc)
  }

  # Species covs
  spc <- x@speciesCovs
  if(!is.null(spc)){
    # M x S covs
    spc_site <- sapply(spc, function(x) identical(dim(x), c(M, S))) 
    spc_site <- spc[spc_site]
    if(length(spc_site) > 0){
      spc_site <- lapply(spc_site, function(x){
        x[i,,drop=FALSE]
      })
    }
    # M x J x S covs
    spc_obs <- sapply(spc, function(x) identical(dim(x), c(M, J, S)))
    spc_obs <- spc[spc_obs]
    if(length(spc_obs) > 0){
      spc_obs <- lapply(spc_obs, function(x){
        x[i,,,drop=FALSE]
      })
    }
    new_spc <- c(spc_site, spc_obs)
  } else {
    new_spc <- NULL
  }

  umf <- x
  umf@y <- ylist[[1]]
  umf@ylist <- ylist
  umf@siteCovs <- siteCovs
  umf@obsCovs <- obsCovs
  umf@speciesCovs <- new_spc
  umf
})

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
  J <- obsNum(object@data)
  new_object <- object
  newform <- multispeciesFormula(object@formula)
  new_object@formula <- newform$formula
  new_object@data <- process_multispecies_umf(object@data, newform$covs)
  new_object <- as(new_object, "unmarkedFitOccu")

  pr <- predict(new_object, type=type, newdata=newdata, backTransform=backTransform,
                na.rm=na.rm, appendData=appendData, level=level, re.form=re.form, ...)

  if(type == "state"){
    inds <- split(1:nrow(pr), rep(1:length(object@data@ylist), each=M))
  } else if(type == "det"){
    inds <- split(1:nrow(pr), rep(1:length(object@data@ylist), each=M*J))
  }
  names(inds) <- names(object@data@ylist)
  lapply(inds, function(x){ 
         out <- pr[x,,drop=FALSE]
         rownames(out) <- NULL
         out
  })
})
