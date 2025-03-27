# Functions for the book Applied Hierarchical Modeling in Ecology (AHM)
# Marc Kery & Andy Royle, Academic Press, 2016.

# pcount.spHDS - section 9.8.4 p540 --- SEE ERRATA

# Function fits spatial hierarchical distance sampling model
# unmarked model fitting function introduced in Section 9.8.4)
# update 4/28/2017

pcount.spHDS<- function (formula, data, K, mixture = c("P", "NB"), 
  starts = NULL, method = "BFGS", se = TRUE, ...){

  # Check arguments
  forms <- split_formula(formula)
  check_no_support(forms)
  if (!is(data, "unmarkedFramePCount")){
    stop("Data is not an unmarkedFramePCount object.")
  }
  mixture <- match.arg(mixture)

  # Build submodels
  fam <- switch(mixture, P = "poisson", NB = "negative_binomial")
  submodels <- unmarkedSubmodelList(
    state = unmarkedSubmodelState(name = "Abundance", short_name = "lam", 
                                  type = "state", formula = forms[[2]], data = data, 
                                  family = fam, link = "log"),

    det = unmarkedSubmodelDet(name = "Detection", short_name = "p", 
                              type = "det", formula = forms[[1]], data = data, 
                              family = "binomial", link = "logit")
  )

  if(mixture == "NB"){
    submodels['alpha'] <- unmarkedSubmodelScalar(name = "Dispersion", 
                                                 short_name = "alpha", 
                                                 type = "alpha", link = "log")
  }

  # Build response object
  response <- unmarkedResponseCount(data, submodels, Kmax = K)

  inputs <- nll_inputs(response, submodels, engine = "R")

  # Fit model
  fit <- fit_model(nll_pcount.spHDS_R, inputs = inputs, submodels = submodels,
                   starts = starts, method = method, se = se, ...)

  new("unmarkedFitPCount", fitType = "pcount", call = match.call(),
      formula = formula, data = data, sitesRemoved = removed_sites(response),
      estimates = fit$submodels, AIC = fit$AIC, opt = fit$opt, 
      negLogLike = fit$opt$value, nllFun = fit$nll, 
      K = response@Kmax, mixture = mixture)
}


nll_pcount.spHDS_R <- function(params, inputs){
  with(inputs, {

  beta_state <- params[idx_state[1]:idx_state[2]]
  beta_det <- params[idx_det[1]:idx_det[2]]
  beta_det[1] <- -1 * exp(beta_det[1])

  M <- nrow(y)
  J <- ncol(y)

  k <- 0:Kmax
  k.ik <- rep(k, M)
  k.ijk <- rep(k, M * J)
  y.ij <- as.numeric(t(y))
  y.ijk <- rep(y.ij, each = Kmax + 1)
  ijk <- expand.grid(k = 0:Kmax, j = 1:J, i = 1:M)
  ijk.to.ikj <- with(ijk, order(i, k, j))

  theta.i <- exp(X_state %*% beta_state + offset_state)
  p.ij <- 2 * plogis(X_det %*% beta_det + offset_det)
  
  theta.ik <- rep(theta.i, each = Kmax + 1)
  p.ijk <- rep(p.ij, each = Kmax + 1)
  bin.ijk <- dbinom(y.ijk, k.ijk, p.ijk)

  bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (Kmax + 1), J, byrow = TRUE)
  g.ik <- rowProds(bin.ik.mat)
  if (family_state == 1) { # Poisson
    f.ik <- dpois(k.ik, theta.ik)
  } else if (family_state == 2) { # negbin
    alpha <- exp(params[idx_alpha[1]])
    f.ik <- dnbinom(k.ik, mu = theta.ik, size = alpha)
  }
  dens.i.mat <- matrix(f.ik * g.ik, M, Kmax + 1, byrow = TRUE)
  dens.i <- rowSums(dens.i.mat)
  -sum(log(dens.i), na.rm=TRUE)
  })
}
