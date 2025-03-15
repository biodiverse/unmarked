setClass("unmarkedResponse",
  slots = c(
    y = "ANY",
    Kmin = "numeric"
  ), 
  prototype = list(
    y = numeric(0),
    Kmin = numeric(0)
  )
)

unmarkedResponse <- function(data){
  response <- new("unmarkedResponse", y = data@y)
  #response <- add_missing(response, submodels)
  Kmin <- apply(data@y, 1, function(x){
    ifelse(all(is.na(x)), NA, max(x, na.rm=TRUE))
  })
  response@Kmin <- Kmin
  response
}

#setGeneric("add_missing", function(resp, submod, ...){
#  standardGeneric("add_missing")
#})

#setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelState"),
#  function(resp, submod, ...){
#  
#  resp@y[find_missing(submod),] <- NA
#  resp
#})

#setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelDet"),
#  function(resp, submod, ...){
#  yt <- t(resp@y)
#  yt[find_missing(submod)] <- NA
#  resp@y <- t(yt)
#  resp
#})

#setMethod("add_missing", c("unmarkedResponse", "unmarkedSubmodelList"),
#  function(resp, submod, ...){  
#  for (i in 1:length(submod@submodels)){
#    resp <- add_missing(resp, submod@submodels[[i]])
#  }
#  resp
#})

setMethod("get_TMB_data", "unmarkedResponse", function(object){
  list(y = object@y, Kmin = object@Kmin)
})
