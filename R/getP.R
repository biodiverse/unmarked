#setGeneric("getP", function(object, ...) standardGeneric("getP"))

# Exported method
setMethod("getP", "unmarkedFit", function(object, ...){
  getP_internal(object)
})

# Internal method
setGeneric("getP_internal", function(object) standardGeneric("getP_internal"))

setMethod("getP_internal", "unmarkedFit", function(object){
  M <- numSites(object@data)
  J <- obsNum(object@data)
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- matrix(p, M, J, byrow = TRUE)
  p
})

setMethod("getP_internal", "unmarkedFitColExt", function(object){
  data <- object@data
  M <- numSites(data)
  nY <- data@numPrimary
  J <- obsNum(data)/nY
  p <- predict(object, type="det", level=NULL, na.rm=FALSE)$Predicted
  p <- array(p, c(J, nY, M))
  p <- aperm(p, c(3, 1, 2))
  p <- matrix(p, nrow=M)
  p
})

setMethod("getP_internal", "unmarkedFitDS", function(object){
  cp <- get_dist_prob(object)
  cp
})

setMethod("getP_internal", "unmarkedFitMPois", function(object){
  p <- methods::callNextMethod(object)
  piFun <- object@data@piFun
  pi <- do.call(piFun, list(p = p))
  pi
})

setMethod("getP_internal", "unmarkedFitOccuMS", function(object){
  J <- ncol(object@data@y)
  N <- nrow(object@data@y)
  pred <- predict(object, type='det', level=NULL)
  lapply(pred, function(x) matrix(x$Predicted, nrow=N, ncol=J, byrow=T))
})

setMethod("getP_internal", "unmarkedFitOccuMulti", function(object){
  ylist <- object@data@ylist
  S <- length(ylist)
  N <- nrow(ylist[[1]])
  maxOrder <- object@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
  dm <- getDesign(object@data,object@detformulas,object@stateformulas, maxOrder=maxOrder)
  pred <- predict(object,'det', level=NULL)
  dets <- do.call(cbind,lapply(pred,`[`,,1))
  #ugly mess
  out <- list()
  for (i in 1:S){
    pmat <- array(NA,dim(ylist[[1]]))
    for (j in 1:N){
      ps <- dets[dm$yStart[j]:dm$yStop[j],i]
      not_na <- !is.na(ylist[[i]][j,])
      pmat[j,not_na] <- ps
    }
    out[[i]] <- pmat
  }
  names(out) <- names(ylist)
  out
})

setMethod("getP_internal", "unmarkedFitOccuTTD", function(object){
  N <- nrow(object@data@y)
  lam <- predict(object, 'det', na.rm=FALSE)$Predicted
  tmax <- as.numeric(t(object@data@surveyLength))
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")

  not_na <- !is.na(lam)
  est_p <- rep(NA, length(lam))
  if(tdist == "weibull"){
    k <- exp(coef(object)['k(k)'])
    est_p[not_na] <- stats::pweibull(tmax[not_na], k, 1/lam[not_na])
  } else {
    est_p[not_na] <- stats::pexp(tmax[not_na], lam[not_na])
  }

  matrix(est_p, nrow=N, byrow=TRUE)
})
