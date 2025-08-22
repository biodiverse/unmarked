# This file starts with All_ so it is loaded first by R

# Class unions-----------------------------------------------------------------
setClassUnion("optionalDataFrame", c("data.frame","NULL"))
setClassUnion("optionalMatrix", c("matrix","NULL"))
setClassUnion("optionalNumeric", c("numeric","NULL"))
setClassUnion("optionalCharacter", c("character","NULL"))
setClassUnion("optionalList", c("list","NULL"))
setClassUnion("numericOrLogical", c("numeric", "logical"))
setClassUnion("matrixOrVector", c("matrix","numeric"))

# unmarkedFrame classes--------------------------------------------------------

validunmarkedFrame <- function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors,
               "siteCovData does not have same size number of sites as y.")
    if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
    if(length(errors) == 0)
        TRUE
    else
        errors
}

setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)

## a class for multi-season data

setClass("unmarkedMultFrame",
    representation(numPrimary = "numeric",
        #data frame in site-major, year-minor order describing siteCovs
        yearlySiteCovs = "optionalDataFrame"),
    contains="unmarkedFrame")

## a class for distance sampling data
setClass("unmarkedFrameDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrame",
    validity = function(object) {
        errors <- character(0)
        J <- numY(object)
        db <- object@dist.breaks
        if(J != length(db) - 1)
            errors <- c(errors, "ncol(y) must equal length(dist.breaks)-1")
        if(db[1] != 0)
            errors <- c(errors, "dist.breaks[1] must equal 0")
        if(!is.null(obsCovs(object)))
            "obsCovs cannot be used with distsamp"
        if(length(errors) == 0) TRUE
        else errors
        })

setClassUnion("unmarkedFrameOrNULL", members=c("unmarkedFrame", "NULL"))

setClass("unmarkedFrameOccu",
		contains = "unmarkedFrame")

setClass("unmarkedFrameOccuComm",
         representation(ylist = "list", speciesCovs="optionalList"),
         contains = "unmarkedFrame")

setClass("unmarkedFrameOccuFP",
         representation(
           type = "numeric"),
         contains = "unmarkedFrame")

#Multi-state occupancy
setClass('unmarkedFrameOccuMS',
         representation(
            numStates = "numeric",
            phiOrder = "list"),
         contains = "unmarkedMultFrame")

#Time-to-detection occupancy
setClass("unmarkedFrameOccuTTD",
         representation(
            surveyLength = "matrix"),
          contains = "unmarkedMultFrame")

setClass("unmarkedFrameOccuMulti",
    representation(ylist = "list", fDesign = "matrix"),
    contains = "unmarkedFrame",
    validity = function(object) {
      errors <- character(0)
      M <- nrow(object@y)
      J <- ncol(object@y)
      Ms <- sapply(object@ylist,nrow)
      Js <- sapply(object@ylist,ncol)
      if(length(unique(Ms)) != 1)
        errors <- c(errors, "All species must have same number of sites")
      if(length(unique(Js)) != 1)
        errors <- c(errors, "All species must have same number of observations")
      if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors,
               "siteCovData does not have same size number of sites as y.")
      if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
      if(length(errors) == 0)
        TRUE
      else
        errors
      })

setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")


setClass("unmarkedFrameMPois",
		representation(
			samplingMethod = "character",
			piFun = "character"),
		contains = "unmarkedFrame")


setClass("unmarkedFrameG3",
         contains = "unmarkedMultFrame")

setClass("unmarkedFrameGMM",
    representation(
        piFun = "character",
        samplingMethod = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGPC",
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGOccu", contains = "unmarkedFrameG3")

setClass("unmarkedFrameDailMadsen",
         representation(primaryPeriod = "matrix"),
         contains = "unmarkedMultFrame")

setClass("unmarkedFramePCO",
         contains = "unmarkedFrameDailMadsen")

setClass("unmarkedFrameMMO",
    representation(
        piFun = "character",
        samplingMethod = "character"),
    contains = "unmarkedFrameDailMadsen")

setClass("unmarkedFrameDSO",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrameDailMadsen")

setClass("unmarkedFrameGDR",
  representation(
    yDistance = "matrix",
    yRemoval = "matrix",
    survey = "character",
    dist.breaks = "numeric",
    unitsIn = "character",
    period.lengths = "numeric"
  ),
  contains="unmarkedMultFrame"
)

setClass("unmarkedFrameOccuCOP",
  representation(L = "matrix"),
  contains = "unmarkedFrame",
  validity = function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    y_integers = sapply(object@y, check.integer, na.ignore = T)
    if (!all(y_integers)) {
      errors <- c(errors,
                  paste(
                    "Count detection should be integers. Non-integer values:",
                    paste(object@y[which(!y_integers)], collapse = ', ')
                  ))
    }
    if (!all(all(dim(object@L) == dim(object@y)))){
      errors <- c( errors, paste(
        "L should be a matrix of the same dimension as y, with M =", M,
        "rows (sites) and J =", J, "columns (sampling occasions)."
      ))}
    if (length(errors) == 0) TRUE
    else errors
  }
)

# unmarkedFit classes----------------------------------------------------------

# Class to store actual parameter estimates
setClass("unmarkedEstimate",
    representation(
        name = "character",
		short.name = "character",
        estimates = "numeric",
        covMat = "matrix",
        fixed = "numeric",
        covMatBS = "optionalMatrix",
        invlink = "character",
        invlinkGrad = "character",
        randomVarInfo= "list"),
    validity = function(object) {
        errors <- character(0)
        if(nrow(object@covMat) != length(object@estimates)) {
        errors <- c(errors,
            "Size of covMat does not match length of estimates.")
        }
    if(length(errors) > 0)
        errors
    else
        TRUE
    })

setClass("unmarkedEstimateList",
    representation(estimates = "list"),
    validity = function(object) {
        errors <- character(0)
        for(est in object@estimates) {
            if(!is(est, "unmarkedEstimate")) {
                errors <- c("At least one element of unmarkedEstimateList is not an unmarkedEstimate.")
                break
                }
            }
        if(length(errors) == 0) {
            return(TRUE)
        } else {
            return(errors)
            }
    })

setClass("unmarkedFit",
   representation(fitType = "character",
        call = "call",
        formula = "formula",
        data = "unmarkedFrame",
        sitesRemoved = "numeric",  # vector of indices of removed sites
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
        opt = "list",
        negLogLike = "numeric",
        nllFun = "function",
        bootstrapSamples = "optionalList",
        covMatBS = "optionalMatrix", # list of bootstrap sample fits
        TMB = "optionalList")) #TMB output object

# ---------------------------- CHILD CLASSES ----------------------------

setClass("unmarkedFitDS",
    representation(
        keyfun = "character",
        unitsOut = "character",
        output = "character"),
        contains = "unmarkedFit")

setClass("unmarkedFitGOccu",
    representation(
        formlist = "list"),
    contains = "unmarkedFit")

setClass("unmarkedFitPCount",
    representation(
        K = "numeric",
        mixture = "character"),
    contains = "unmarkedFit")

# This class is not used directly, just used as a base for for PCO, MMO, DSO
setClass("unmarkedFitDailMadsen",
        representation(
            K = "numeric",
            mixture = "character",
            formlist = "list",
            dynamics = "character",
            immigration = "logical",
            fix = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitPCO", contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitMMO", contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitDSO",
        representation(
            keyfun = "character",
            unitsOut = "character",
            output = "character"),
        contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitIDS",
    representation(
        formlist = "list",
        keyfun = "character",
        K = "numeric",
        dataPC = "unmarkedFrameOrNULL",
        dataOC = "unmarkedFrameOrNULL",
        maxDist = "list",
        surveyDurations = "list",
        unitsOut = "character"),
        contains = "unmarkedFit")

setClass("unmarkedFitOccu",
    representation(knownOcc = "logical"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuComm", contains="unmarkedFitOccu")

setClass("unmarkedFitOccuCOP",
         representation(removed_obs = "matrix",
                        formlist = "list"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuPEN",
    representation(
	knownOcc = "logical",
	pen.type = "character",
	lambda = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuPEN_CV",
    representation(
	knownOcc = "logical",
	pen.type = "character",
	lambdaVec = "numeric",
	k = "numeric",
	foldAssignments = "numeric",
	lambdaScores = "numeric",
	chosenLambda = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuFP",
         representation(knownOcc = "logical",
            detformula = "formula",
            FPformula = "formula",
            Bformula = "formula",
            stateformula = "formula",
            type = "numeric"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuMulti",
         representation(
            detformulas = "character",
            stateformulas = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuMS",
         representation(
            detformulas = "character",
            psiformulas = "character",
            phiformulas = "character",
            parameterization = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuTTD",
    representation(
        psiformula = "formula",
        gamformula = "formula",
        epsformula = "formula",
        detformula = "formula"),
    contains = "unmarkedFit")

setClass("unmarkedFitNmixTTD",
         representation(
           stateformula = "formula",
           detformula = "formula",
           K = "numeric"),
         contains = "unmarkedFit")

setClass("unmarkedFitMPois",
    contains = "unmarkedFit")


setClass("unmarkedFitOccuRN",
    representation(
      K = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitColExt",
    representation(
        phi = "matrix",
        psiformula = "formula",
        gamformula = "formula",
        epsformula = "formula",
        detformula = "formula",
        projected = "array",
        projected.mean = "matrix",
        smoothed = "array",
        smoothed.mean = "matrix",
        projected.mean.bsse = "optionalMatrix",
        smoothed.mean.bsse = "optionalMatrix"),
    contains = "unmarkedFit")


setClass("unmarkedFitGMM",
    representation(
        formlist = "list",
        mixture = "character",
        K = "numeric"),
    contains = "unmarkedFit")


setClass("unmarkedFitGDS",
    representation(
        keyfun = "character",
        unitsOut = "character",
        output = "character"),
    contains = "unmarkedFitGMM")

setClass("unmarkedFitGDR", contains = "unmarkedFitGDS")

setClass("unmarkedFitGPC",
    contains = "unmarkedFitGMM")


# Other output objects---------------------------------------------------------
# Linear combinations object
setClass("unmarkedLinComb",
         representation(parentEstimate = "unmarkedEstimate",
                        estimate = "numeric",
                        covMat = "matrix",
                        covMatBS = "optionalMatrix",
                        coefficients = "matrix"))

setClass("unmarkedBackTrans",
         representation(parentLinComb = "unmarkedLinComb",
                        estimate = "numeric",
                        covMat = "matrix",
                        covMatBS = "optionalMatrix"))

setClassUnion("linCombOrBackTrans", c("unmarkedLinComb", "unmarkedBackTrans"))

setClass("profile", representation(prof = "matrix"))

setClass("unmarkedRanef",
    representation(post = "array"))

setClass("unmarkedFitList",
    representation(fits = "list"),
    validity = function(object) {
        fl <- object@fits
        umf1 <- getData(fl[[1]])
        y1 <- fl_getY(fl[[1]])
        dataTest <- sapply(fl, function(x) isTRUE(all.equal(umf1, getData(x))))
        yTest <- sapply(fl, function(x) isTRUE(all.equal(y1, fl_getY(x))))
        if(!all(dataTest)) {
            stop("Data are not the same among models. Make sure you use the same unmarkedFrame object for all models.")
            }
        else if(!all(yTest)) {
            stop("Data are not the same among models due to missing covariate values. Consider removing NAs before analysis.")
            }
        TRUE
        }
    )

setClass("unmarkedModSel",
    representation(
        Full = "data.frame",
        Names = "matrix"
        )
    )

setClass("parboot",
         representation(call = "call",
                        t0 = "numeric",
                        t.star = "matrix"))

setClass("unmarkedCrossVal",
    representation(stats = "data.frame",
                   summary = "data.frame",
                   method = "character",
                   folds = "numeric",
                   holdoutPct = "numeric"),
    validity=function(object){
      errors <- character(0)
      hp <- object@holdoutPct
      if(hp<0|hp>1){
        errors <- c(errors,"holdoutPct must be between 0 and 1")
      }
    }
)

setClass("unmarkedCrossValList",
    representation(stats_list="list",
                   method = "character",
                   folds="numeric",
                   holdoutPct="numeric",
                   sort="character")
)

setClass("unmarkedPower",
  representation(call="call", data="unmarkedFrame", M="numeric",
                 J="numeric", T="numeric", coefs="list", estimates="list",
                 alpha="numeric", nulls="list")
)

setClass("unmarkedPowerList", representation(powerAnalyses="list"))

setClass("unmarkedPostSamples",
         representation(numSites="numeric",
                        numPrimary="numeric",
                        nsims="numeric",
                        samples="array")
         )
