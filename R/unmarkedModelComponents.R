setClass("unmarkedModelComponents",
  slots = c(
    response = "unmarkedResponse",
    submodels = "unmarkedSubmodelList",
    auxiliary = "list"
  )
)

unmarkedModelComponents <- function(response, submodels, auxiliary = list()){
  
  new("unmarkedModelComponents", response = response,
      submodels = submodels, auxiliary = auxiliary)
}

setMethod("get_TMB_data", "unmarkedModelComponents", function(object){
  c(get_TMB_data(object@response),
    get_TMB_data(object@submodels),
    get_TMB_data(object@auxiliary))
})

setGeneric("submodelList", function(object) standardGeneric("submodelList"))

setMethod("submodelList", "unmarkedModelComponents", function(object){
  object@submodels
})

setGeneric("submodelList<-", function(object, value){
  standardGeneric("submodelList<-")
})

setMethod("submodelList<-", "unmarkedModelComponents", function(object, value){
  object@submodels <- value
  object
})
