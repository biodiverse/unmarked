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

# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula, data, sitesRemoved,
    estimates, AIC, opt, negLogLike, nllFun)
{
    umfit <- new("unmarkedFit", fitType = fitType, call = call,
        formula = formula, data = data, sitesRemoved = sitesRemoved,
        estimates = estimates, AIC = AIC, opt = opt,
        negLogLike = negLogLike,
        nllFun = nllFun)
    return(umfit)
}

# ---------------------------- CHILD CLASSES ----------------------------


setClass("unmarkedFitDS",
    representation(
        keyfun = "character",
        unitsOut = "character",
        output = "character"),
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

setClass("unmarkedFitOccu",
    representation(knownOcc = "logical"),
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


setClass("unmarkedFitGPC",
    contains = "unmarkedFitGMM")



# -------------------------- Show and Summary ----------------------------


setMethod("show", "unmarkedFit", function(object)
{
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    show(object@estimates)
    cat("AIC:", object@AIC,"\n")
    if(!identical(object@opt$convergence, 0L))
        warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
})




setMethod("summary", "unmarkedFit", function(object)
{
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    summaryOut <- summary(object@estimates)
    cat("AIC:", object@AIC,"\n")
    cat("Number of sites:", sampleSize(object))
    if(length(object@sitesRemoved) > 0)
        cat("\nID of sites removed due to NA:", object@sitesRemoved)
    cat("\noptim convergence code:", object@opt$convergence)
    cat("\noptim iterations:", object@opt$counts[1], "\n")
    if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
    cat("Bootstrap iterations:", length(object@bootstrapSamples), "\n\n")
    invisible(summaryOut)
})



setMethod("summary", "unmarkedFitDS", function(object)
{
    out <- callNextMethod()
    cat("Survey design: ", object@data@survey, "-transect", sep="")
    cat("\nDetection function:", object@keyfun)
    cat("\nUnitsIn:", object@data@unitsIn)
    cat("\nUnitsOut:", object@unitsOut, "\n\n")
    invisible(out)
})




# Compute linear combinations of estimates in unmarkedFit objects.
setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "matrixOrVector"),
    function(obj, coefficients, type, offset = NULL, re.form=NULL)
{
    stopifnot(!missing(type))
    stopifnot(type %in% names(obj))
    estimate <- obj@estimates[type]
    linearComb(estimate, coefficients, offset, re.form)
})


setMethod("backTransform", "unmarkedFit", function(obj, type)
{
    est <- obj[type]
    if(length(est@estimates) == 1) {
        lc <- linearComb(est, 1)
        return(backTransform(lc))
    } else {
        stop('Cannot directly backTransform an unmarkedEstimate with length > 1.')
        }
})


setMethod("[", "unmarkedFit",
          function(x, i, j, drop) {
              x@estimates[i]
          })


setMethod("names", "unmarkedFit",
          function(x) {
              names(x@estimates)
          })



# ---------------------- coef, vcov, and SE ------------------------------


setMethod("coef", "unmarkedFit",
    function(object, type, altNames = TRUE, fixedOnly=TRUE)
{
    if(missing(type)) {
        co <- lapply(object@estimates@estimates,
            function(x) coef(x, altNames=altNames, fixedOnly=fixedOnly))
        names(co) <- NULL
        co <- unlist(co)
    } else {
        co <- coef(object[type], altNames=altNames, fixedOnly=fixedOnly)
        }
    co
})


setMethod("vcov", "unmarkedFit",
    function (object, type, altNames = TRUE, method = "hessian", fixedOnly=TRUE, ...)
{
    method <- match.arg(method, c("hessian", "nonparboot"))
    if(.hasSlot(object, "TMB") && !is.null(object@TMB)) method <- "TMB"
    switch(method,
           hessian = {
            if (is.null(object@opt$hessian)) {
                stop("Hessian was not computed for this model.")
            }
            v <- solve(hessian(object))
        },
        nonparboot = {
            if (is.null(object@bootstrapSamples)) {
                stop("No bootstrap samples have been drawn. Use nonparboot first.")
                }
            v <- object@covMatBS
        },
        TMB = {
          return(vcov_TMB(object, type, fixedOnly))
        })
    rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
    if (missing(type)) {
        return (v)
    } else {
      object[type]
        inds <- .estimateInds(object)[[type]]
        return (v[inds, inds, drop = FALSE])
        }
})

## A helper function to return a list of indices for each estimate type
##
.estimateInds <- function(umf) {
  ## get length of each estimate
  estimateLengths <- sapply(umf@estimates@estimates, function(est) {
    length(coef(est))
  })
  ## recurse function to generate list of indices
  estimateInds <- function(type) {
    if(type==1) {
      return(list(seq(length=estimateLengths[1])))
    } else {
      prev.list <- estimateInds(type-1)
      prev.max <- max(prev.list[[type-1]])
      return(c(prev.list, list(seq(prev.max+1, prev.max +
                                   estimateLengths[type]))))
    }
  }
  retlist <- estimateInds(length(estimateLengths))
  names(retlist) <- names(umf)
  retlist
}

setMethod("vcov", "unmarkedFitOccuMulti",
    function (object, type, altNames = TRUE, method = "hessian", ...)
{

  pen <- object@call[["penalty"]]
  if(is.null(pen)) pen <- 0
  if(pen>0) method <- "nonparboot"
  callNextMethod(object, type, altNames, method=method, ...)

})


setMethod("SE", "unmarkedFit", function(obj,...)
{
    v <- vcov(obj,...)
    sqrt(diag(v))
})



setMethod("logLik", "unmarkedFit", function(object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    ll <- -object@negLogLike
    #attr(ll, "df") <- length(coef(object))
    #class(ll) <- "logLik"
    return(ll)
})



setMethod("LRT", c(m1="unmarkedFit", m2="unmarkedFit"), function(m1, m2)
{
    ll1 <- logLik(m1)
    ll2 <- logLik(m2)
    chisq <- 2 * abs(ll1 - ll2)
    DF <- abs(length(coef(m1)) - length(coef(m2)))
    pval <- pchisq(chisq, DF, lower.tail=FALSE)
    return(data.frame(Chisq=chisq, DF = DF, 'Pr(>Chisq)' = pval,
        check.names=F))
})





setMethod("confint", "unmarkedFit", function(object, parm, level = 0.95,
    type, method = c("normal", "profile"))
{
    method <- match.arg(method)
    if(missing(type))
        stop(paste("Must specify type as one of (", paste(names(object), collapse=", "),").",sep=""))
    if(missing(parm))
        parm <- object[type]@fixed
    if(method == "normal") {
        callGeneric(object[type],parm = parm, level = level)
    } else {
        nllFun <- nllFun(object)
        ests <- mle(object)
        nP <- length(parm)
        ci <- matrix(NA, nP, 2)

        ## create table to match parm numbers with est/parm numbers
        types <- names(object)
        numbertable <- data.frame(type = character(0), num = numeric(0))
        for(i in seq(length=length(types))) {
            length.est <- length(object[i]@estimates)
            numbertable <- rbind(numbertable, data.frame(type =
                rep(types[i], length.est), num = seq(length=length.est)))
            }
        parm.fullnums <- which(numbertable$type == type &
            numbertable$num %in% parm)

        for(i in seq(length=nP)) {
            cat("Profiling parameter",i,"of",nP,"...")
            se <- SE(object[type])
            whichPar <- parm.fullnums[i]
            ci[i,] <- profileCI(nllFun, whichPar=whichPar, MLE=ests,
                interval=ests[whichPar] + 10*se[i]*c(-1,1), level=level)
            cat(" done.\n")
            }
        rownames(ci) <- names(coef(object[type]))[parm]
        colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
        return(ci)
        }
})


setMethod("profile", "unmarkedFit",
    function(fitted, type, parm, seq)
{
    stopifnot(length(parm) == 1)
    MLE <- mle(fitted)
    SE(fitted[type])
    nPar <- length(mle(fitted))
    nll <- nllFun(fitted)

    ## create table to match parm numbers with est/parm numbers
    types <- names(fitted)
    numbertable <- data.frame(type = character(0), num = numeric(0))
    for(i in seq(length=length(types))) {
        length.est <- length(fitted[i]@estimates)
        numbertable <- rbind(numbertable,
        data.frame(type = rep(types[i], length.est),
            num = seq(length=length.est)))
        }
    parm.fullnums <- which(numbertable$type == type &
        numbertable$num == parm)

    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, parm.fullnums, value)
        mleRestricted <- optim(rep(0,nPar), fixedNLL)$value
        mleRestricted
        }
        prof.out <- sapply(seq, f)
        prof.out <- cbind(seq, prof.out)
        new("profile", prof = prof.out)
})



setMethod("hessian", "unmarkedFit",
    function(object)
{
    object@opt$hessian
})


setMethod("update", "unmarkedFit", function(object, ..., evaluate=TRUE){
  inp <- list(...)
  cl <- object@call
  for (i in 1:length(inp)){
    # Issue here with creating a new list element that is NULL
    # Need to do it in two steps
    cl[[names(inp)[i]]] <- NA # first create empty element
    cl[[names(inp)[i]]] <- inp[[i]] # then fill it
  }

  if(!evaluate) return(cl)

  eval(cl, parent.frame(2))
})


setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))
setMethod("sampleSize", "unmarkedFit", function(object) {
    data <- getData(object)
    M <- numSites(data)
    M <- M - length(object@sitesRemoved)
    M
})


setGeneric("getData", function(object) standardGeneric("getData"))
setMethod("getData", "unmarkedFit",function(object) {
    object@data
})


setGeneric("nllFun", function(object) standardGeneric("nllFun"))
setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)

setGeneric("mle", function(object) standardGeneric("mle"))
setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setClass("profile", representation(prof = "matrix"))

setGeneric("smoothed",
    function(object, mean=TRUE) standardGeneric("smoothed"))
setMethod("smoothed","unmarkedFitColExt",
    function(object, mean) {
        if(mean) object@smoothed.mean
        else object@smoothed
    })

setGeneric("projected",
    function(object, mean=TRUE) standardGeneric("projected"))
setMethod("projected","unmarkedFitColExt", function(object, mean) {
    if(mean) object@projected.mean
    else object@projected
})

setMethod("plot", c("profile", "missing"), function(x) {
    plot(x@prof[,1], x@prof[,2], type = "l")
})


setMethod("plot", c(x = "unmarkedFit", y = "missing"), function(x, y, ...)
{
    r <- residuals(x)
    e <- fitted(x, na.rm = FALSE)
    plot(e, r, ylab = "Residuals", xlab = "Predicted values")
    abline(h = 0, lty = 3, col = "gray")
})

setMethod("plot", c(x = "unmarkedFitOccuMulti", y = "missing"), function(x, y, ...)
{
  r <- do.call(rbind,residuals(x))
  e <- do.call(rbind,fitted(x))
  plot(e, r, ylab = "Residuals", xlab = "Predicted values")
  abline(h = 0, lty = 3, col = "gray")
})


setMethod("hist", "unmarkedFitDS", function(x, lwd=1, lty=1, ...) {
    ymat <- getY(getData(x))
    dbreaks <- getData(x)@dist.breaks
    nb <- length(dbreaks)
    mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
    distances <- rep(mids, times=colSums(ymat))
    h <- hist(distances, plot=F, breaks=dbreaks)
    key <- x@keyfun
    survey <- x@data@survey
    switch(key,
           halfnorm = {
               sigma <- exp(coef(x, type="det"))
               if(length(sigma) > 1)
               stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- 2 * integrate(dnorm, dbreaks[1],
                                               dbreaks[nb], sd=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) 2 * dnorm(x, mean=0, sd=sigma),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhn, dbreaks[1], dbreaks[nb],
                                           sigma=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhn(r, sigma=sigma),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           exp = {		# This doesn't work on example fm4
               rate <- exp(coef(x, type="det"))
               if(length(rate) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxexp(x, rate=rate),
                               min(dbreaks),
                               max(dbreaks), add=T, lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drexp(r, rate=rate),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      })
           },
           hazard = {
               shape <- exp(coef(x, type="det"))
               scale <- exp(coef(x, type="scale"))
               if(length(shape) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxhaz(x, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhaz(r, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           uniform = {
               switch(survey,
                      line = {
                          plot(h, freq=F, ...)
                          abline(h=1/max(dbreaks), lwd=lwd, lty=lty)
                      },
                      point = {
                          plot(h, freq=F, ...)
                          plot(function(r) (pi*r^2) / (pi*max(dbreaks)^2),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      }
                      )}
           )
})




# ----------------------- CHILD CLASS METHODS ---------------------------

# Extract detection probs
setGeneric("getP", function(object, ...) standardGeneric("getP"))
setGeneric("getFP", function(object, ...) standardGeneric("getFP"))
setGeneric("getB", function(object, ...) standardGeneric("getB"))


setMethod("getFP", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  formula <- object@formula
  detformula <- object@detformula
  stateformula <- object@stateformula
  FPformula <- object@FPformula
  Bformula <- object@Bformula
  umf <- object@data
  designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
  type = object@type
  y <- designMats$y
  U <- designMats$U
  U.offset <- designMats$U.offset
  if (is.null(U.offset))
    U.offset <- rep(0, nrow(U))
  M <- nrow(y)
  J <- ncol(y)
  fpars <- coef(object, type = "fp")
  f <- plogis(U %*% fpars + U.offset)
  f <- matrix(f, M, J, byrow = TRUE)
  if (type[1]!=0){
    f[,1:type[1]] = 0
  }
  return(f)
})

setMethod("getB", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  formula <- object@formula
  detformula <- object@detformula
  stateformula <- object@stateformula
  FPformula <- object@FPformula
  Bformula <- object@Bformula
  umf <- object@data
  designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
  y <- designMats$y
  W <- designMats$W
  W.offset <- designMats$W.offset
  if (is.null(W.offset))
    W.offset <- rep(0, nrow(W))
  M <- nrow(y)
  J <- ncol(y)
  type = object@type
  if (type[3]!=0){
    bpars <- coef(object, type = "b")
  b <- plogis(W %*% bpars + W.offset)
  b <- matrix(b, M, J, byrow = TRUE)
  }
  if (type[3]==0){
    b <- matrix(0, M, J)
  }
  return(b)
})






# Should this return p or pi. Right now it's pi without phi.
setMethod("getP", "unmarkedFitGDS",
    function(object, na.rm = TRUE)
{
#    browser()
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    Xdet <- designMats$Xdet
    Xdet.offset <- designMats$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    ppars <- coef(object, type = "det")
    db <- umf@dist.breaks
    w <- diff(db)
    survey <- umf@survey
    key <- object@keyfun
    tlength <- umf@tlength

    u <- a <- matrix(NA, M, J)
    switch(survey,
    line = {
        for(i in 1:M) {
            a[i,] <- tlength[i] * w
            u[i,] <- a[i,] / sum(a[i,])
            }
        },
    point = {
        for(i in 1:M) {
            a[i, 1] <- pi*db[2]^2
            for(j in 2:J)
                a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
            u[i,] <- a[i,] / sum(a[i,])
            }
        })


    cp <- array(NA, c(M, J, T))
    switch(key,
        halfnorm = {
            sigma <- exp(Xdet %*% ppars + Xdet.offset)
            sigma <- matrix(sigma, M, T, byrow=TRUE)
            for(i in 1:M) {
                for(t in 1:T) {
                    if(is.na(sigma[i,t])) next
                    switch(survey,
                    line = {
                        f.0 <- 2 * dnorm(0, 0, sd=sigma[i, t])
                        int <- 2 * (pnorm(db[-1], 0, sd=sigma[i, t]) -
                            pnorm(db[-(J+1)], 0, sd=sigma[i, t]))
                        cp[i,,t] <- int / f.0 / w
                        },
                    point = {
                        for(j in 1:J) {
                            cp[i, j, t] <- integrate(grhn, db[j], db[j+1],
                                sigma=sigma[i, t], rel.tol=1e-4)$value *
                                2 * pi / a[i, j]
                            }
                        })
                    cp[i,,t] <- cp[i,,t] * u[i,]
                    }
                }
            },
        exp = {
            rate <- exp(Xdet %*% ppars + Xdet.offset)
            rate <- matrix(rate, M, T, byrow=TRUE)
            for(i in 1:M) {
                for(t in 1:T) {
                if(is.na(rate[i,t])) next
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(gxexp, db[j], db[j+1],
                            rate=rate[i,t], rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(grexp, db[j], db[j+1],
                            rate=rate[i,t], rel.tol=1e-4)$value *
                            2 * pi / a[i, j]
                        }
                    })
                cp[i,,t] <- cp[i,,t] * u[i,]
                }}
            },
        hazard = {
            shape <- exp(Xdet %*% ppars + Xdet.offset)
            shape <- matrix(shape, M, T, byrow=TRUE)
            scale <- exp(coef(object, type="scale"))
            for(i in 1:M) {
                for(t in 1:T) {
                if(is.na(shape[i,t])|is.na(scale)) next
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(gxhaz, db[j], db[j+1],
                            shape=shape[i,t], scale=scale,
                            rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(grhaz, db[j], db[j+1],
                            shape = shape[i,t], scale=scale,
                            rel.tol=1e-4)$value * 2 * pi / a[i, j]
                    }})
                cp[i,,t] <- cp[i,,t] * u[i,]
                }}
            },
	uniform = {
#            browser()
            cp[] <- u
        })
    cp <- matrix(cp, nrow=M)
    return(cp)
})


setMethod("getP", "unmarkedFitDSO",
    function(object, na.rm = TRUE)
{
    umf <- getData(object)
    y <- getY(umf)
    M <- numSites(umf)

    keep <- 1:M
    if(na.rm & length(object@sitesRemoved)>0){
      keep <- keep[-object@sitesRemoved]
    }

    M <- length(keep)
    T <- umf@numPrimary
    J <- ncol(y) / T

    sig <- matrix(NA, M, T)
    if(object@keyfun != "uniform"){
      sig <- predict(object, type="det", level = NULL, na.rm=na.rm)$Predicted
      sig <- matrix(sig, M, T, byrow=TRUE)
    }

    scale <- 0.0
    if(object@keyfun == "hazard"){
      scale <- backTransform(object, type="scale")@estimate
    }

    db <- umf@dist.breaks
    w <- diff(db)
    ua <- getUA(umf)
    u <- ua$u[keep,]; a <- ua$a[keep,]

    cp <- array(NA, c(M, J, T))
    for (i in 1:M){
      for (t in 1:T){
        cp[i,,t] <- getDistCP(object@keyfun, sig[i,t], scale, umf@survey,
                              db, w, a[i,], u[i,])
      }
    }
    matrix(cp, nrow=M)
})


setMethod("getP", "unmarkedFitMMO", function(object, na.rm = TRUE)
{

  umf <- object@data

  D <- getDesign(umf, object@formula, na.rm=na.rm)
  beta <- coef(object, type='det')
  off <- D$Xp.offset
  if(is.null(off)) off <- rep(0, nrow(D$Xp))
  plong <- plogis(D$Xp %*% beta + off)

  M <- nrow(D$y)
  T <- umf@numPrimary
  J <- ncol(getY(umf)) / T

  pmat <- aperm(array(plong, c(J,T,M)), c(3,1,2))

  pout <- array(NA, c(M,J,T))
  for (t in 1:T){
    pout[,,t] <- do.call(umf@piFun, list(p=pmat[,,t]))
  }
  matrix(aperm(pout,c(2,3,1)), M, J*T, byrow=TRUE)

})


setMethod("getP", "unmarkedFitPCO", function(object, na.rm = TRUE)
{
    umf <- object@data
    D <- getDesign(umf, object@formula, na.rm = na.rm)
    y <- D$y
    Xp <- D$Xp
    Xp.offset <- D$Xp.offset
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*J)
    ppars <- coef(object, type = "det")
    p <- plogis(Xp %*% ppars + Xp.offset)
    p <- matrix(p, M, J*T, byrow = TRUE)
    return(p)
})



setMethod("getP", "unmarkedFitGMM",
    function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- object@formlist$pformula
    piFun <- object@data@piFun
    umf <- object@data
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xdet <- D$Xdet
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))

    M <- nrow(y)
    T <- object@data@numPrimary
    R <- numY(object@data) / T
    J <- obsNum(object@data) / T

    ppars <- coef(object, type = "det")
    p <- plogis(Xdet %*% ppars + Xdet.offset)
    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, J, T))
    p <- aperm(p, c(1,3,2))

    cp <- array(as.numeric(NA), c(M, T, R))
    for(t in 1:T) cp[,t,] <- do.call(piFun, list(matrix(p[,t,], M, J)))
    cp <- aperm(cp, c(1,3,2))
    cp <- matrix(cp, nrow=M, ncol=numY(object@data))

    return(cp)
})


setMethod("getP", "unmarkedFitGPC",
    function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- object@formlist$pformula
    umf <- object@data
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xdet <- D$Xdet
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))
    R <- nrow(y)
    T <- object@data@numPrimary
    J <- ncol(y) / T
    ppars <- coef(object, type = "det")
    p <- plogis(Xdet %*% ppars + Xdet.offset)
    p <- matrix(p, nrow=R, byrow=TRUE)
    return(p)
})

#Y extractors for unmarkedFit objects
setMethod("getY", "unmarkedFit", function(object) object@data@y)
setMethod("getY", "unmarkedFitOccu", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitOccuRN", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitColExt", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitOccuMulti", function(object) {
            object@data@ylist
})
