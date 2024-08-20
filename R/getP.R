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

setMethod("getP_internal", "unmarkedFitMPois", function(object){
  p <- methods::callNextMethod(object)
  piFun <- object@data@piFun
  pi <- do.call(piFun, list(p = p))
  pi
})
