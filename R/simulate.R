replace_estimates <- function(object, new_ests){
  for (i in 1:length(new_ests)){
    est <- object@estimates@estimates[[names(new_ests)[i]]]@estimates
    stopifnot(length(est) == length(new_ests[[i]]))
    object@estimates@estimates[[names(new_ests)[i]]]@estimates <- new_ests[[i]]
  }
  object
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
          stop("Random effect covariates must be specified as factors with guide argument", call.=FALSE)
        }
        b <- stats::rnorm(length(levels(lvldata)), 0, old_coefs[signame])
        names(b) <- rep(paste0("b_",i), length(b))
        new_coefs <- c(new_coefs, b)
        coefs[[i]] <- new_coefs
      }
    }
  }
  coefs
}

