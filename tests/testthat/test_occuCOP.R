context("occuCOP fitting function")
skip_on_cran()

COPsimul <- function(psi = 0.5,
                         lambda = 1,
                         M = 100,
                         J = 5) {
  
  z_i <- sample(
    x = c(0, 1),
    size = M,
    prob = c(1 - psi, psi),
    replace = T
  )
  
  y = matrix(rpois(n = M * J, lambda = lambda), nrow = M, ncol = J) * z_i
  
  return(y)
}


# COPsimulfit <- function(psi = 0.5,
#                         lambda = 1,
#                         M = 100,
#                         J = 5,
#                         obsLength = NULL,
#                         ...) {
#   occuCOP(data = unmarkedFrameCOP(
#     y = simulDataCOP(
#       psi = psi,
#       lambda = lambda,
#       M = M,
#       J = J
#     ),
#     obsLength = obsLength
#   ), ...)
# }

test_that("unmarkedFrameCOP is constructed correctly", {
  set.seed(123)
  
  # Simulate data
  M = 100
  J = 5
  y = COPsimul(psi = .5,
                   lambda = 1,
                   M = M,
                   J = J)
  obsLength = y * 0 + 1
  
  psiCovs = data.frame(
    "psiNum" = rnorm(M),
    "psiInt" = as.integer(rpois(n = M, lambda = 5)),
    "psiBool" = sample(c(T, F), size = M, replace = T),
    "psiChar" = sample(letters[1:5], size = M, replace = T),
    "psiFactUnord" = factor(sample(
      letters[5:10], size = M, replace = T
    )),
    "psiFactOrd" = sample(factor(c("A", "B", "C"), ordered = T), size =
                            M, replace = T)
  )
  
  lambdaCovs = list(
    "lambdaNum" = matrix(
      rnorm(M * J), 
      nrow = M, ncol = J
    ),
    "lambdaInt" = matrix(
      as.integer(rpois(n = M * J, lambda = 1e5)),
      nrow = M, ncol = J
    ),
    "lambdaBool" = matrix(
      sample(c(T, F), size = M * J, replace = T),
      nrow = M, ncol = J
    ),
    "lambdaChar" = matrix(
      sample(letters[1:5], size = M * J, replace = T),
      nrow = M, ncol = J
    ),
    "lambdaFactUnord" = matrix(
      factor(sample(letters[5:10], size = M * J, replace = T)),
      nrow = M, ncol = J
    ),
    "lambdaFactOrd" = matrix(
      sample(factor(c("A", "B", "C"), ordered = T), size = M * J, replace = T),
      nrow = M, ncol = J
    )
  )

  
  # Creating a unmarkedFrameCOP object
  expect_warning(umf <- unmarkedFrameCOP(y = y))
  expect_no_error(umf <- unmarkedFrameCOP(y = y, obsLength = obsLength))
  
  
  # Create subsets
  expect_no_error(umf_sub_i <- umf[1:3, ])
  expect_no_error(umf_sub_j <- umf[, 1:2])
  expect_no_error(umf_sub_ij <- umf[1:3, 1:2])
  
  # unmarkedFrameCOP organisation ----------------------------------------------
  expect_true(inherits(umf, "unmarkedFrameCOP"))
  expect_equivalent(numSites(umf_sub_i), 3)
  expect_equivalent(obsNum(umf_sub_j), 2)
  expect_equivalent(numSites(umf_sub_ij), 3)
  expect_equivalent(obsNum(umf_sub_ij), 2)
  
  # unmarkedFrameCOP display ---------------------------------------------------
  
  # print method
  expect_no_error(print(umf))
  expect_no_warning(print(umf))
  expect_no_message(print(umf))
  expect_no_error(print(umf_sub_i))
  expect_no_error(print(umf_sub_j))
  expect_no_error(print(umf_sub_ij))
  expect_no_error(print(umf[1,]))
  expect_no_error(print(umf[,1]))
  expect_no_error(print(umf[1,1]))
  
  # summary method for unmarkedFrameCOP
  expect_no_error(summary(umf))
  expect_no_error(summary(umf_sub_ij))
  
  # plot method for unmarkedFrameCOP
  expect_no_error(plot(umf))
  expect_no_error(plot(umf_sub_ij))
  
  
  # Input handling: covariates -------------------------------------------------
  expect_no_error(umfCovs <- unmarkedFrameCOP(
    y = y,
    obsLength = obsLength,
    siteCovs = psiCovs,
    obsCovs = lambdaCovs
  ))
  expect_no_error(print(umfCovs))
  expect_no_error(summary(umfCovs))
  expect_no_error(plot(umfCovs))
  
  # Input handling: NA ---------------------------------------------------------
  
  # NA should not pose problems when creating the unmarkedFrameCOP object.
  # The warnings and potential errors will be displayed when fitting the model.
  # Except when y only contains NA: then there's an error.
  
  ## NA in y
  yNA <- y
  yNA[1:2,] <- NA
  yNA[3:5, 3:4] <- NA
  yNA[,3] <- NA
  expect_no_error(umfNA <- unmarkedFrameCOP(y = yNA, obsLength = obsLength))
  expect_no_error(print(umfNA))
  expect_no_error(summary(umfNA))
  expect_no_error(plot(umfNA))
  
  ## NA in obsLength
  obsLengthNA <- obsLength
  obsLengthNA[1:2,] <- NA
  obsLengthNA[3:5, 3:4] <- NA
  obsLengthNA[,5] <- NA
  expect_no_error(umfNA <- unmarkedFrameCOP(y = y, obsLength = obsLengthNA))
  expect_no_error(print(umfNA))
  expect_no_error(summary(umfNA))
  expect_no_error(plot(umfNA))
  
  ## NA also in covariates
  psiCovsNA <- psiCovs
  psiCovsNA[1:5,] <- NA
  psiCovsNA[c(8,10,12), 3] <- NA
  psiCovsNA[,1] <- NA
  lambdaCovsNA <- lambdaCovs
  lambdaCovsNA[[1]][1:5,] <- NA
  lambdaCovsNA[[1]][,3] <- NA
  lambdaCovsNA[[3]][,5] <- NA
  expect_no_error(umfCovsNA <- unmarkedFrameCOP(
    y = yNA,
    obsLength = obsLengthNA,
    siteCovs = psiCovsNA,
    obsCovs = lambdaCovsNA
  ))
  expect_no_error(print(umfCovsNA))
  expect_no_error(summary(umfCovsNA))
  expect_no_error(plot(umfCovsNA))
  
  ## all NA in y
  yallNA <- y
  yallNA[1:M, 1:J] <- NA
  expect_error(unmarkedFrameCOP(y = yallNA, obsLength = obsLength))
  
  # Input handling: error and warnings -----------------------------------------
  # Error when there are decimals in y
  y_with_decimals = y
  y_with_decimals[1, 1] = .5
  expect_error(unmarkedFrameCOP(y = y_with_decimals, obsLength = obsLength))
  
  # Warning if y is detection/non-detection instead of count
  y_detec_nodetec = (y > 0) * 1
  expect_warning(unmarkedFrameCOP(y = y_detec_nodetec, obsLength = obsLength))
  
  # Error if the dimension of obsLength is different than that of y
  expect_error(unmarkedFrameCOP(y = y, obsLength = obsLength[1:2, 1:2]))
})


test_that("occuCOP can fit simple models", {
  # Simulate data
  set.seed(123)
  M = 100
  J = 5
  y = COPsimul(psi = .5,
                   lambda = 1,
                   M = M,
                   J = J)
  obsLength = y * 0 + 1
  
  # Create umf
  umf <- unmarkedFrameCOP(y = y, obsLength = obsLength)
  
  # Fitting options ----
  
  # With default parameters
  expect_message(fit_default <- occuCOP(data = umf))
  expect_no_error(occuCOP(data = umf, psiformula = ~ 1, lambdaformula = ~ 1, psistarts = 0, lambdastarts = 0))
  
  # With chosen starting points
  expect_no_error(occuCOP(data = umf, 
                          psiformula = ~ 1, lambdaformula = ~ 1, 
                          psistarts = qlogis(.7),
                          lambdastarts = log(0.1)))
  expect_error(occuCOP(data = umf, 
                       psiformula = ~ 1, lambdaformula = ~ 1, 
                       psistarts = qlogis(c(0.7, 0.5), lambdastarts = 0)))
  expect_error(occuCOP(data = umf, 
                       psiformula = ~ 1, lambdaformula = ~ 1, 
                       lambdastarts = log(c(1, 2)), psistarts = 0))
  
  # With different options
  expect_no_error(occuCOP(data = umf, method = "Nelder-Mead", psistarts = 0, lambdastarts = 0))
  expect_error(occuCOP(data = umf, method = "ABC", psistarts = 0, lambdastarts = 0))
  
  expect_no_error(occuCOP(data = umf, se = F, psistarts = 0, lambdastarts = 0))
  expect_error(occuCOP(data = umf, se = "ABC"))
  
  expect_no_error(occuCOP(data = umf, engine = "R", psistarts = 0, lambdastarts = 0))
  expect_error(occuCOP(data = umf, engine = "julia", psistarts = 0, lambdastarts = 0))
  
  expect_no_error(occuCOP(data = umf, na.rm = F, psistarts = 0, lambdastarts = 0))
  expect_error(occuCOP(data = umf, na.rm = "no", psistarts = 0, lambdastarts = 0))
  
  # Looking at at COP model outputs ----
  expect_is(fit_default, "unmarkedFitCOP")
  
  ## backTransform
  expect_no_error(backTransform(fit_default, type = "psi"))
  expect_no_error(backTransform(fit_default, type = "lambda"))
  expect_error(backTransform(fit_default, type = "state"))
  expect_error(backTransform(fit_default, type = "det"))
  expect_is(backTransform(fit_default, type = "psi"), "unmarkedBackTrans")
  
  ## predict with newdata = fit@data
  expect_no_error(predict(object = fit_default, type = "psi"))
  expect_no_error(predict(object = fit_default, type = "lambda"))
  expect_error(predict(object = fit_default, type = "state"))
  
  ## predict with newdata = 1
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame(1),
    type = "psi"
  ))
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame(1),
    type = "lambda"
  ))
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame("a"=1:5,"b"=10:14),
    type = "psi"
  ))
  
  # Fitting accurately ----
  
  ## psi = 0.50, lambda = 1 ----
  psi_test = .5
  lambda_test = 1
  fit_accur <- occuCOP(data = unmarkedFrameCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    obsLength = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate, 
    psi_test, 
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate, 
    lambda_test, 
    tol = 0.05
  )
  
  ## psi = 0.25, lambda = 5 ----
  psi_test = 0.25
  lambda_test = 5
  fit_accur <- occuCOP(data = unmarkedFrameCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    obsLength = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate, 
    psi_test, 
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate, 
    lambda_test, 
    tol = 0.05
  )
  
  ## psi = 0.75, lambda = 0.5 ----
  psi_test = 0.75
  lambda_test = 0.5
  fit_accur <- occuCOP(data = unmarkedFrameCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    obsLength = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate, 
    psi_test, 
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate, 
    lambda_test, 
    tol = 0.05
  )

  # With NAs ----
  yNA <- y
  yNA[1,] <- NA
  yNA[3, 1] <- NA
  yNA[4, 3] <- NA
  yNA[, 5] <- NA
  expect_no_error(umfNA <- unmarkedFrameCOP(y = yNA, obsLength = obsLength))
  
  expect_warning(fit_NA <- occuCOP(data = umfNA, psistarts = 0, lambdastarts = 0))
  expect_error(occuCOP(data = umfNA, psistarts = 0, lambdastarts = 0, na.rm = F))
  
  
  #TODO
})


test_that("occuCOP can fit models with covariates", {
  
})


test_that("We can simulate COP data", {
  
  # From scratch ----
  
  # With no covariates
  expect_no_error(simulate(
    "COP",
    formulas = list(psi =  ~ 1, lambda =  ~ 1),
    coefs = list(
      psi = c(intercept = 0),
      lambda = c(intercept = 0)
    ),
    design = list(M = 100, J = 100)
  ))
  
  # With quantitative covariates
  expect_no_error(simulate(
    "COP",
    formulas = list(psi =  ~ elev, lambda =  ~ rain),
    coefs = list(
      psi = c(intercept = qlogis(.5), elev = -0.5),
      lambda = c(intercept = log(3), rain = -1)
    ),
    design = list(M = 100, J = 5)
  ))
  
  # With guides
  expect_no_error(simulate(
    "COP",
    formulas = list(psi =  ~ elev, lambda =  ~ rain),
    coefs = list(
      psi = c(intercept = qlogis(.5), elev = -0.5),
      lambda = c(intercept = log(3), rain = -1)
    ),
    design = list(M = 100, J = 5),
    guide = list(elev=list(dist=rnorm, mean=12, sd=0.5))
  ))
  
  # With qualitative covariates
  expect_no_error(simulate(
    "COP",
    formulas = list(psi =  ~ elev + habitat, lambda =  ~ 1),
    coefs = list(
      psi = c(
        intercept = qlogis(.2),
        elev = -0.5,
        habitatB = .5,
        habitatC = .8
      ),
      lambda = c(intercept = log(3))
    ),
    design = list(M = 100, J = 5),
    guide = list(habitat = factor(levels = c("A", "B", "C")))
  ))
  
  # From unmarkedFitCOP ----
  #TODO
})

