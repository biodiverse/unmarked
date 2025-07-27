process_site_index <- function(umf, i){
  if(is.logical(i)){
    stopifnot(length(i) == numSites(umf))
    keep <- 1:numSites(umf)
    keep <- keep[i]
  } else if(is.numeric(i)){
    stopifnot(all(abs(i) %in% 1:numSites(umf)))
    stopifnot(all(i>0) | all(i<0))
    if(all(i<0)){
      keep <- 1:numSites(umf)
      keep <- keep[!keep %in% abs(i)]
    } else {
      keep <- i
    }
  }
  keep
}

process_obs_index <- function(umf, j){
  if(methods::.hasSlot(umf, "numPrimary")){
    if(umf@numPrimary > 1){
      stop("Can't subset observations in multi-period umf", call.=FALSE)
    }
  }
  if(is.logical(j)){
    stopifnot(length(j) == obsNum(umf))
    keep <- 1:obsNum(umf)
    keep <- keep[j]
  } else if(is.numeric(j)){
    stopifnot(all(abs(j) %in% 1:obsNum(umf)))
    stopifnot(all(j>0) | all(j<0))
    if(all(j<0)){
      keep <- 1:obsNum(umf)
      keep <- keep[!keep %in% abs(j)]
    } else {
      keep <- j
    }
  }
  keep
}

process_period_index <- function(umf, j){
  stopifnot(methods::.hasSlot(umf, "numPrimary"))
  T <- umf@numPrimary
  stopifnot(T > 1)
  if(is.logical(j)){
    stopifnot(length(j) == T)
    keep <- 1:T
    keep <- keep[j]
  } else if(is.numeric(j)){
    stopifnot(all(abs(j) %in% 1:T))
    stopifnot(all(j>0) | all(j<0))
    if(all(j < 0)){
      keep <- 1:T
      keep <- keep[!keep %in% abs(j)]
    } else {
      keep <- j
    }
  }
  keep
}

subset_site <- function(umf, i){
 
  M <- numSites(umf)
  J <- obsNum(umf)

  # Process selection
  keep <- process_site_index(umf, i)

  # Subset y
  y <- getY(umf)[keep,,drop=FALSE]
  umf@y <- y

  # Subset site covs
  sc <- siteCovs(umf)
  if(!is.null(sc)){
    siteCovs(umf) <- sc[keep,,drop=FALSE]
  }

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_obs <- rep(1:M, each = J)
    oc_new <- lapply(keep, function(x){
      oc[keep_obs == x,,drop=FALSE]
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset yearly site covs if they exist
  if(methods::.hasSlot(umf, "yearlySiteCovs")){
    ysc <- yearlySiteCovs(umf)
    keep_per <- rep(1:M, each = umf@numPrimary)
    if(!is.null(ysc)){
      ysc_new <- lapply(keep, function(x){
        ysc[keep_per == x,,drop=FALSE]
      })
      yearlySiteCovs(umf) <- do.call(rbind, ysc_new)
    }
  }

  umf
}

subset_period <- function(umf, j){
  if(obsNum(umf) != numY(umf)){
    stop("Can't subset multinomial observations", call.=FALSE)
  }
  M <- numSites(umf)
  T_old <- umf@numPrimary
  J <- numY(umf) / T_old
  R <- obsNum(umf) / T_old

  # Process selection
  keep <- process_period_index(umf, j)

  umf@numPrimary <- length(keep)

  # Subset y
  keep_per <- rep(1:T_old, each = J)
  y <- lapply(keep, function(x){
    getY(umf)[,keep_per == x,drop=FALSE]
  })
  umf@y <- do.call(cbind, y)

  # Subset yearly site covs
  ysc <- yearlySiteCovs(umf)
  if(!is.null(ysc)){
    keep_site <- rep(1:M, each = T_old)
    ysc_new <- lapply(1:M, function(m){
      ysc_sub <- ysc[keep_site == m,,drop=FALSE]
      ysc_sub <- ysc_sub[keep,,drop=FALSE]
    })
    yearlySiteCovs(umf) <- do.call(rbind, ysc_new)
  }

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_site <- rep(1:M, each = obsNum(umf))
    keep_per <- rep(1:T_old, each = R)

    oc_new <- lapply(1:M, function(m){
      oc_sub <- oc[keep_site == m,,drop=FALSE]
      oc_sub <- lapply(keep, function(x){
        oc_sub[keep_per == x,,drop=FALSE] 
      })
      do.call(rbind, oc_sub)
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset obsToY
  # Only works with non-multinomial
  obsToY(umf) <- diag(numY(umf)) 

  umf
}

subset_obs <- function(umf, j){
  if(obsNum(umf) != numY(umf)){
    stop("Can't subset multinomial observations", call.=FALSE)
  }
  M <- numSites(umf)
  J <- obsNum(umf)

  # Process selection
  keep <- process_obs_index(umf, j)

  # Subset y
  y <- getY(umf)[,keep,drop=FALSE]
  umf@y <- y

  # Subset obs covs
  oc <- obsCovs(umf)
  if(!is.null(oc)){
    keep_site <- rep(1:M, each = J)
    oc_new <- lapply(1:M, function(m){
      oc_sub <- oc[keep_site == m,,drop=FALSE]
      oc_sub <- oc_sub[keep,,drop=FALSE]
    })
    obsCovs(umf) <- do.call(rbind, oc_new)
  }

  # Subset obsToY
  # This only works with non-multinomial
  obsToY(umf) <- diag(length(keep))

  umf
}

setClassUnion("numericOrLogical", c("numeric", "logical"))

setMethod("[", c("unmarkedFrame", "numericOrLogical", "missing", "missing"),
  function(x, i){
  subset_site(x, i)
})

setMethod("[", c("unmarkedFrame", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  subset_obs(x, j)
})

setMethod("[", c("unmarkedMultFrame", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  subset_period(x, j)
})

setMethod("[", c("unmarkedFrame", "numericOrLogical", "numericOrLogical", "missing"),
  function(x, i, j){
  x <- x[i,]
  x[,j]
})

# unmarkedFrame-specific site-selector methods---------------------------------

setMethod("[", c("unmarkedFrameDS", "numericOrLogical", "missing", "missing"),
  function(x, i){
  idx <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  if(x@survey == "line"){
    x@tlength <- x@tlength[idx] 
  }
  x
})

setMethod("[", c("unmarkedFrameGDS", "numericOrLogical", "missing", "missing"),
  function(x, i){
  idx <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  if(x@survey == "line"){
    x@tlength <- x@tlength[idx] 
  }
  x
})

setMethod("[", c("unmarkedFrameDailMadsen", "numericOrLogical", "missing", "missing"),
  function(x, i){
  idx <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  x@primaryPeriod <- x@primaryPeriod[idx,,drop=FALSE]
  x
})

setMethod("[", c("unmarkedFrameDSO", "numericOrLogical", "missing", "missing"),
  function(x, i){
  idx <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  if(x@survey == "line"){
    x@tlength <- x@tlength[idx] 
  }
  x
})

setMethod("[", c("unmarkedFrameOccuTTD", "numericOrLogical", "missing", "missing"),
  function(x, i){
  idx <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  x@surveyLength <- x@surveyLength[idx,,drop=FALSE]
  x
})

setMethod("[", c("unmarkedFrameOccuMulti", "numericOrLogical", "missing", "missing"),
  function(x, i){
  keep <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  x@ylist <- lapply(x@ylist, function(s) s[keep,,drop=FALSE])
  x
})

setMethod("[", c("unmarkedFrameOccuMulti", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  keep <- process_obs_index(x, j)
  x <- methods::callNextMethod(x, i, j)
  x@ylist <- lapply(x@ylist, function(s) s[,keep,drop=FALSE])
  x
})

setMethod("[", c("unmarkedFrameOccuCOP", "numericOrLogical", "missing", "missing"),
  function(x, i){
  keep <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  x@L <- x@L[keep,,drop=FALSE]
  x
})

setMethod("[", c("unmarkedFrameOccuCOP", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  keep <- process_obs_index(x, j)
  x <- methods::callNextMethod(x, i, j)
  x@L <- x@L[,keep,drop=FALSE]
  x
})

setMethod("[", c("unmarkedFrameGDR", "numericOrLogical", "missing", "missing"),
  function(x, i){
  keep <- process_site_index(x, i)
  x <- methods::callNextMethod(x, i)
  x@yDistance <- x@yDistance[keep,,drop=FALSE]
  x@yRemoval <- x@yRemoval[keep,,drop=FALSE]
  x
})

setMethod("[", c("unmarkedFrameGDR", "missing", "numericOrLogical", "missing"),
  function(x, i, j){
  keep <- process_period_index(x, j)
  T <- x@numPrimary
  Jdist <- ncol(x@yDistance) / T
  Jrem <- ncol(x@yRemoval) / T
  x <- methods::callNextMethod(x, i, j)

  # Subset y matrices
  keep_per <- rep(1:T, each = Jdist)
  yDistance <- lapply(keep, function(k){
    x@yDistance[,keep_per == k,drop=FALSE]
  })
  x@yDistance <- do.call(cbind, yDistance)

  keep_per <- rep(1:T, each = Jrem)
  yRemoval <- lapply(keep, function(k){
    x@yRemoval[,keep_per == k,drop=FALSE]
  })
  x@yRemoval <- do.call(cbind, yRemoval)
  x
})

setMethod("head", "unmarkedFrame", function(x, n) {
  if(missing(n)) n <- 6
  x[1:n,]
})
