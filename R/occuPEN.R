occuPEN <- function(formula, data, knownOcc = numeric(0), starts = NULL,
                    method = "BFGS", engine = c("C", "R"),
		                lambda = 0, pen.type = c("Bayes","Ridge","MPLE"), ...){

  engine <- match.arg(engine)
  pen.type <- match.arg(pen.type)
  out <- occu(formula = formula, data = data, knownOcc = knownOcc, starts = starts,
              method = method, engine = engine, se = FALSE, lambda = lambda,
              pen.type = pen.type, ...)
  out@call <- match.call()
  out
}


occuPEN_CV <- function(formula, data, knownOcc = numeric(0), starts = NULL,
                        method = "BFGS", engine = c("C", "R"),
		                    lambdaVec = c(0,2^seq(-4,4)), pen.type = c("Bayes", "Ridge"),
		                    k = 5, foldAssignments = NA,...){

  if(!is(data, "unmarkedFrameOccu"))
        stop("Data is not an unmarkedFrameOccu object.")

  pen.type <- match.arg(pen.type)

  if (length(lambdaVec)==1) stop("Must provide more than one lambda for cross-validation.")

  engine <- match.arg(engine, c("C", "R"))

  # Fit 1 iteration just to create the submodels, handle missing data etc.
  template_fit <- occu(formula = formula, data = data, knownOcc = knownOcc,
                       starts = starts, method = "SANN", engine = engine,
                       control = list(maxit=1), se = FALSE)

  response <- unmarkedResponseBinary(data, template_fit@estimates)
  has_obs <- response@Kmin == 1
  
  M <- numSites(data)
  J <- obsNum(data)

  # user-supplied foldAssignments
  if(!(length(foldAssignments)==1 & is.na(foldAssignments)[1])) {
    if(k != length(unique(foldAssignments))){
      stop("Value of k does not match number of folds indicated in foldAssignments.")
    }
  } else { # create foldAssignments
    # attempt to include sites with and without observations in each fold
    foldAssignments = c(1:M)
    idxsWithObs = which(has_obs)
    idxsWoObs = which(!has_obs)
    if (length(idxsWithObs)>0 & length(idxsWoObs)>0) {
      foldAssignments[idxsWithObs] <-
        sample(rep(1:k,ceiling(length(idxsWithObs)/k))[1:length(idxsWithObs)])
      foldAssignments[idxsWoObs] <-
        sample(rep(1:k,ceiling(length(idxsWoObs)/k))[1:length(idxsWoObs)])
    } else if (k <= M) {
      foldAssignments = sample(rep(1:k,ceiling(M/k)))[1:M]
    } else {
      stop("k>M. More folds than sites creates folds. Specify a smaller k.")
    }
  }
  foldNames <- unique(foldAssignments)

  lambdaScores <- rep(0, length(lambdaVec)) # score by held-out likelihood

  for (f in 1:k) {
    fold <- foldNames[f]
    occuTrain <- data[which(foldAssignments!=fold),] # train on NOT this fold
    occuTest <- data[which(foldAssignments==fold),] # test on this fold

    for (la in 1:length(lambdaVec)) {
      # First fit the model with training data to get estimates
      fit_train <- occuPEN(formula, occuTrain, starts = starts,
                           lambda = lambdaVec[la], pen.type = pen.type)
      ests_train <- coef(fit_train)

      # Create the nll function for the test data by fitting
      # one iteration (for speed)
      fit_test_nll <- occu(formula, occuTest, starts = starts, se = FALSE,
                           method = "SANN", control = list(maxit = 1))@nllFun

      # Apply the estimates to the test data
      lambdaScores[la] <- lambdaScores[la] + fit_test_nll(ests_train)
    }
  }

  bestLambda <- lambdaVec[which.min(lambdaScores)]

  final <- occuPEN(formula, data, knownOcc = knownOcc, starts=starts, 
                   engine = engine, lambda=bestLambda, pen.type=pen.type, ...)

  known_occ <- rep(FALSE, numSites(data))
  known_occ[knownOcc] <- TRUE

  new("unmarkedFitOccuPEN_CV", fitType = "occu", call = match.call(),
      formula = formula, data = data, sitesRemoved = removed_sites(response),
      estimates = final@estimates, AIC = final@AIC, opt = final@opt,
      negLogLike = final@negLogLike, nllFun = final@nllFun, knownOcc = known_occ,
		  pen.type = pen.type, lambdaVec = lambdaVec, k = k, 
      foldAssignments = foldAssignments, lambdaScores = lambdaScores, 
      chosenLambda = bestLambda)
}
