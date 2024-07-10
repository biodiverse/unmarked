setMethod("simulate", "unmarkedFrame",
  function(object, nsim = 1, seed = NULL, model = NULL, coefs = NULL,
           quiet = FALSE, ...){
  object <- y_to_zeros(object)
  fit <- get_fit(object, model, ...)
  coefs <- check_coefs(coefs, fit, quiet = quiet)
  coefs <- generate_random_effects(coefs, fit)
  fit <- replace_estimates(fit, coefs)
  sims <- simulate(fit, nsim)
  lapply(sims, function(x) replaceY(object, x))
})

setGeneric("y_to_zeros", function(object, ...){
  standardGeneric("y_to_zeros")
})

# Other fit-specific methods at the bottom of the file
setMethod("y_to_zeros", "unmarkedFrame", function(object, ...){
  object@y[] <- 0
  object
})

get_fit <- function(object, model, ...){
  fun <- get_fitting_function(object, model)
  fun(..., data = object, method = "SANN",
      control=list(maxit=0), se=FALSE)
}

setGeneric("get_fitting_function", function(object, model, ...){
  standardGeneric("get_fitting_function")
})

# Other fit-specific methods at the bottom of the file
setMethod("get_fitting_function", "unmarkedFrameOccu",
          function(object, model, ...){
 
  if(!(identical(model, occuRN) | identical(model, occu) | identical(model, occuPEN))){
    stop("model argument must be occu, occuRN, or occuPEN", call.=FALSE)
  }
  model
})

check_coefs <- function(coefs, fit, name = "coefs", quiet = FALSE){
  required_subs <- names(fit@estimates@estimates)
  required_coefs <- lapply(fit@estimates@estimates, function(x) names(x@estimates))

  formulas <- sapply(names(fit), function(x) get_formula(fit, x))

  # If there are random effects, adjust the expected coefficient names
  # to remove the b vector and add the grouping covariate name
  rand <- lapply(formulas, lme4::findbars)
  if(!all(sapply(rand, is.null))){
    stopifnot(all(required_subs %in% names(formulas)))
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    if(!all(sapply(rvar, length)<2)){
      stop("Only 1 random effect per parameter is supported", call.=FALSE)
    }
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- required_coefs[[i]]
        new_coefs <- old_coefs[!grepl("b_", old_coefs, fixed=TRUE)]
        new_coefs <- c(new_coefs, signame)
        required_coefs[[i]] <- new_coefs
      }
    }
  }
  required_lens <- lapply(required_coefs, length)

  dummy_coefs <- lapply(required_coefs, function(x){
                    out <- rep(0, length(x))
                    names(out) <- x
                    out
                  })

  if(is.null(coefs)){
    cat(name, "should be a named list of vectors, with the following structure
        (replace 0s with your values):\n\n")
    print(dummy_coefs)
    stop(paste("Specify", name, "argument as shown above"), call.=FALSE)
  }

  for (i in 1:length(required_subs)){
    if(!required_subs[i] %in% names(coefs)){
      stop(paste0("Missing required list element '",
                  required_subs[i], "' in ", name, " list"), call.=FALSE)
    }

    sub_coefs <- coefs[[required_subs[i]]]

    if(!quiet){
      message(paste0("Assumed parameter order for ", required_subs[i], ":\n",
                  paste(required_coefs[[i]], collapse=", ")))
    }
    
    if(length(sub_coefs) != required_lens[i]){
      stop(paste0("Entry '",required_subs[[i]], "' in ", name, " list must be length ",
                  required_lens[[i]]), call.=FALSE)
    }
    names(coefs[[required_subs[i]]]) <- required_coefs[[i]]

  }
  coefs[required_subs]
}

generate_random_effects <- function(coefs, fit){
  required_subs <- names(fit@estimates@estimates)
  formulas <- sapply(names(fit), function(x) get_formula(fit, x))
  rand <- lapply(formulas, lme4::findbars)
  if(!all(sapply(rand, is.null))){
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- coefs[[i]]
        new_coefs <- old_coefs[names(old_coefs)!=signame]

        # Find levels of factor variable
        if(signame %in% names(siteCovs(fit@data))){
          lvldata <- siteCovs(fit@data)[[signame]]
        } else if(signame %in% names(obsCovs(fit@data))){
          lvldata <- obsCovs(fit@data)[[signame]]
        } else if(methods::.hasSlot(fit@data, "yearlySiteCovs") && signame %in% names(yearlySiteCovs(fit@data))){
          lvldata <- yearlySiteCovs(fit@data)[[signame]]
        } else {
          stop("Random effect covariate missing from data", call.=FALSE)
        }

        if(!is.factor(lvldata)){
          stop("Random effect covariates must be specified as factors", call.=FALSE)
        }
        sigma <- old_coefs[signame]
        if(sigma <= 0){
          stop("estimate for random effect represents sigma and must be positive",
               call.=FALSE)
        }
        b <- stats::rnorm(length(levels(lvldata)), 0, sigma)
        names(b) <- rep(paste0("b_",i), length(b))
        new_coefs <- c(new_coefs, b)
        coefs[[i]] <- new_coefs
      }
    }
  }
  coefs
}

replace_estimates <- function(object, new_ests){
  for (i in 1:length(new_ests)){
    est <- object@estimates@estimates[[names(new_ests)[i]]]@estimates
    stopifnot(length(est) == length(new_ests[[i]]))
    object@estimates@estimates[[names(new_ests)[i]]]@estimates <- new_ests[[i]]
  }
  object
}


# y_to_zeros-------------------------------------------------------------------

setMethod("y_to_zeros", "unmarkedFrameOccuMulti", function(object, ...){
  newy <- lapply(object@ylist, function(x){
    x[] <- 0
    x
  })
  object@ylist <- newy
  object
})

setMethod("y_to_zeros", "unmarkedFrameGDR", function(object, ...){
  object@yDistance[] <- 0
  object@yRemoval[] <- 0
  object
})


# get_fitting_function---------------------------------------------------------

setMethod("get_fitting_function", "unmarkedFrameDS", 
          function(object, model, ...){
  if(!missing(model) && identical(model, IDS)) stop("IDS not supported", call.=FALSE)
  distsamp
})

setMethod("get_fitting_function", "unmarkedFrameDSO", 
          function(object, model, ...){
  distsampOpen
})

setMethod("get_fitting_function", "unmarkedFrameGDR", 
          function(object, model, ...){
  gdistremoval
})

setMethod("get_fitting_function", "unmarkedFrameGDS",
          function(object, model, ...){
  gdistsamp
})

setMethod("get_fitting_function", "unmarkedFrameGMM", 
          function(object, model, ...){
  gmultmix
})

setMethod("get_fitting_function", "unmarkedFrameGPC", 
          function(object, model, ...){
  gpcount
})

setMethod("get_fitting_function", "unmarkedFrameGOccu",
          function(object, model, ...){
  goccu
})

setMethod("get_fitting_function", "unmarkedFrameOccuCOP",
          function(object, model, ...){
  occuCOP
})

setMethod("get_fitting_function", "unmarkedFrameOccuFP",
          function(object, model, ...){
  occuFP
})

setMethod("get_fitting_function", "unmarkedFrameOccuMS",
          function(object, model, ...){
  occuMS
})

setMethod("get_fitting_function", "unmarkedFrameOccuTTD",
          function(object, model, ...){
  if(!(identical(model, occuTTD) | identical(model, nmixTTD))){
    stop("model argument must be occuTTD or nmixTTD", call.=FALSE)
  }
  model
})

setMethod("get_fitting_function", "unmarkedFrameOccuMulti",
          function(object, model, ...){
  occuMulti
})

setMethod("get_fitting_function", "unmarkedFrameMPois",
          function(object, model, ...){
  multinomPois
})

setMethod("get_fitting_function", "unmarkedFrameMMO", 
          function(object, model, ...){
  multmixOpen
})

setMethod("get_fitting_function", "unmarkedFramePCount",
          function(object, model, ...){
  pcount
})

setMethod("get_fitting_function", "unmarkedFramePCO", 
          function(object, model, ...){
  pcountOpen
})

setMethod("get_fitting_function", "unmarkedMultFrame", 
          function(object, model, ...){
  colext
})
