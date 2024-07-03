setClass("unmarkedPower",
  representation(call="call", data="unmarkedFrame", M="numeric",
                 J="numeric", T="numeric", coefs="list", estimates="list",
                 alpha="numeric", nulls="list")
)

setGeneric("powerAnalysis", function(object, ...){
  standardGeneric("powerAnalysis")
})

# unmarkedFrame method
# TODO: random effects handling
# TODO: parallel processing
setMethod("powerAnalysis", "unmarkedFrame",
          function(object, model = NULL, effects = NULL, alpha = 0.05,
                   nsim = 100, nulls = NULL, ...){

  test_data <- y_to_zeros(object)
  test_fit <- get_fit(test_data, model, ...)
  effects <- check_coefs(effects, test_fit, name = "effects")

  data_sims <- simulate(object, nsim = nsim, model = model, coefs = effects,
                        quiet = TRUE, ...)

  powerAnalysis_internal(object, model, data_sims, effects, alpha, nulls, ...)
})

# list of unmarkedFrames (pre-simulated) method
setMethod("powerAnalysis", "list",
          function(object, model = NULL, effects = NULL, alpha = 0.05,
                   nsim = length(object), nulls = NULL, ...){
    
  data1 <- object[[1]]
  stopifnot(inherits(data1, "unmarkedFrame"))
  stopifnot(all(sapply(object, function(x) identical(class(data1), class(x)))))
  stopifnot(nsim <= length(object))
  object <- object[1:nsim]

  test_data <- y_to_zeros(data1)
  fit <- get_fit(test_data, model, ...)
  effects <- check_coefs(effects, fit, name = "effects")

  powerAnalysis_internal(data1, model, object, effects, alpha, nulls, ...)
})

powerAnalysis_internal <- function(object, model, data_sims, 
                                   effects, alpha, nulls, ...){
  
  fun <- get_fitting_function(object, model)
  test_fit <- get_fit(data_sims[[1]], model, ...)
  modname <- test_fit@fitType
  modname <- "test"

  if(is.null(nulls)){
    nulls <- effects
    nulls <- lapply(nulls, function(x){
                      x[] <- 0
                      x
                    })
  } else {
    nulls <- check_coefs(nulls, test_fit, name = "nulls")
  }

  sum_dfs <- pbapply::pblapply(data_sims, function(x){
              fit <- fun(..., data = x)
              get_summary_df(fit, effects, nulls)
          })

  sites <- numSites(object)
  primaryPeriods <- ifelse(methods::.hasSlot(object, "numPrimary"),
                           object@numPrimary, 1)
  occasions <- ncol(object@y) / primaryPeriods

  new("unmarkedPower", call=call(modname), data=object, 
      M=sites, J=occasions, T=primaryPeriods, 
      coefs=effects, estimates=sum_dfs, alpha=alpha, nulls=nulls)
}

get_summary_df <- function(fit, effects, nulls){
  n_est <- length(fit@estimates@estimates)
  est_names <- names(fit@estimates@estimates)
  all_est <- lapply(1:n_est, function(i){
    utils::capture.output(out <- summary(fit@estimates@estimates[[i]]))
    out <- out[,1:2]
    out <- cbind(submodel=est_names[i], param=rownames(out), out)
    rownames(out) <- NULL
    out
  })
  all_est <- do.call(rbind, all_est)
  # TODO: Remove random effects
  all_est$Effect <- unlist(effects[est_names])
  all_est$Null <- unlist(nulls[est_names])

  for (i in 1:nrow(all_est)){
    # wald and diff_dir in utils.R
    all_est$P[i] <- wald(all_est$Estimate[i], all_est$SE[i], all_est$Null[i])
    all_est$Direct[i] <- diff_dir(all_est$Estimate[i], all_est$Effect[i],
                                          all_est$Null[i])
  }
  all_est
}

wald <- function(est, se, null_hyp=NULL){
  if(is.null(null_hyp) || is.na(null_hyp)) null_hyp <- 0
  Z <- (est-null_hyp)/se
  2*pnorm(abs(Z), lower.tail = FALSE)
}

diff_dir <- function(est, hyp, null_hyp=NULL){
  if(is.null(null_hyp) || is.na(null_hyp)) null_hyp <- 0
  dif <- est - null_hyp
  dif_hyp <- hyp - null_hyp
  dif * dif_hyp > 0
}

setMethod("summary", "unmarkedPower", 
          function(object, alpha, showIntercepts=FALSE, ...){

  out <- object@estimates[[1]][,c(1,2,5,6)]
  names(out)[1:2] <- c("Submodel", "Parameter")

  if(missing(alpha)){
    alpha <- object@alpha
  }
  stopifnot(alpha >= 0 & alpha <= 1)

  for (i in 1:nrow(out)){
    pcrit <- sapply(object@estimates, function(x) x$P[i]) < alpha
    direct <- sapply(object@estimates, function(x) x$Direct[i])
    ests <- sapply(object@estimates, function(x) x$Estimate[i])

    out$Power[i] <- mean(pcrit & direct, na.rm=TRUE)
    out$`Type S`[i] <- sum(pcrit & !direct, na.rm=TRUE) / sum(pcrit, na.rm=TRUE)
    out$`Type M`[i] <- NA

    # Calculate Type M
    # Adjust for null != 0
    diff_null <- out$Effect[i] - out$Null[i]
    # Don't calculate if effect size is 0
    if(diff_null != 0){
      diffs <- ests - out$Null[i]
      out$`Type M`[i] <- mean(abs(diffs[pcrit]), na.rm=TRUE) / abs(diff_null) 
    }
  }

  if(!showIntercepts){
    no_int <- !grepl("(Intercept)",out$Parameter, fixed=TRUE)
    multi_interact <- grepl("\\[.*:.*\\]", out$Parameter)
    keep <- no_int | multi_interact
    if(sum(keep) > 0) out <- out[keep,,drop=FALSE]
  }

  out
})

setMethod("show", "unmarkedPower", function(object){
  cat("Model:", deparse(object@call[[1]]))
  cat("\nSites:", object@M)
  cat("\nPrimary Periods:", object@T)
  cat("\nOccasions:", object@J)
  cat("\nalpha:", object@alpha)
  cat("\n\n")

  cat("Power Statistics:\n")
  sumtab <- summary(object)
  sumtab$Power <- round(sumtab$Power, 3)
  sumtab$`Type M` <- round(sumtab$`Type M`, 3)
  sumtab$`Type S` <- round(sumtab$`Type S`, 3)

  if(all(sumtab$Null == 0)){
    sumtab$Null <- NULL
  }

  print(sumtab, row.names=FALSE)
})

setMethod("plot", c(x="unmarkedPower", y="missing"),
          function(x, y, alpha, showIntercepts = FALSE, ...){
  if(missing(alpha)) alpha <- x@alpha
  stopifnot(alpha >= 0 & alpha <= 1)
  pars <- x@estimates[[1]]$param
  inds <- 1:length(pars)
  if(!showIntercepts){
    no_int <- !grepl("(Intercept)", pars, fixed=TRUE)
    multi_interact <- grepl("\\[.*:.*\\]", pars)
    keep <- no_int | multi_interact
    if(sum(keep) > 0) inds <- which(keep)
    #inds <- which(pars != "(Intercept)")
  }
  
  if(length(inds) > 1){
    old_ask <- devAskNewPage()
    devAskNewPage(TRUE)
  }
  sapply(inds, function(i) plot_power(x, i, alpha=alpha, ...))
  if(length(inds) > 1) devAskNewPage(old_ask)
  invisible()
})

plot_power <- function(object, ind, alpha, ...){

  submod <- object@estimates[[1]]$submodel[ind]
  param <- object@estimates[[1]]$param[ind]
  parname <- paste(submod, param, sep=" / ")
  effect <- object@estimates[[1]]$Effect[ind]

  ests <- sapply(object@estimates, function(x) x$Estimate[ind])
  pval <- sapply(object@estimates, function(x) x$P[ind])
  direct <- sapply(object@estimates, function(x) x$Direct[ind])
  
  if(missing(alpha)){
    alpha <- object@alpha
  }
  stopifnot(alpha >= 0 & alpha <= 1)

  idx <- 1:length(ests)
  plot(idx, ests, pch=19, col="gray",
       xlab="Simulation", ylab="Estimated effect size", main=parname, ...)
  points(idx[pval < alpha & direct], ests[pval < alpha & direct], pch=19, col='red')
  points(idx[pval < alpha & !direct], ests[pval < alpha & !direct], pch=19, col='blue')
  abline(h = effect, lty=2, lwd=1.3)

  sig_direct <- ests[pval < alpha & direct]
  if(length(sig_direct > 0)){
    abline(h = mean(sig_direct), lty=2, lwd=1.3, col='red')
  }

  legend('bottomright', pch=19, col=c("gray", "red", "blue"), 
         legend=c("Non-significant", "Significant", "Sig & wrong sign"))
  legend('bottomleft', lty=2, col=c("black", "red"),
         legend=c("True effect size", "Avg significant effect"))
  invisible()
}

# unmarkedPowerlist stuff------------------------------------------------------

setGeneric("unmarkedPowerList", function(object, ...){
             standardGeneric("unmarkedPowerList")})
