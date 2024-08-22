context("occuFP fitting function")
skip_on_cran()

test_that("occuFP model can be fit",{
  set.seed(123)
  n = 100
  o = 10
  o1 = 5
  y = matrix(0,n,o)
  p = .7
  r = .5
  fp = 0.05
  y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
  y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
  y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
  type <- c((o-o1),o1,0)
  site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
  occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
  site <- data.frame(habitat = site)
  occ <- list(METH = occ)
  umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

  m1 <- occuFP(detformula = ~ METH, FPformula = ~1,
               stateformula = ~ habitat, data = umf1)
  expect_equal(names(m1), c("state","det","fp"))
  expect_message(fl <- fitList(m1,m1))
  expect_is(fl,"unmarkedFitList")
  expect_equal(length(fl@fits), 2)

  pr <- predict(m1, "fp")
  expect_equal(dim(pr), c(1000, 4))

  # Check update
  new_fit <- update(m1, stateformula=~1, data=umf1[1:10,])
  expect_equal(length(coef(new_fit, "state")), 1)
  expect_equal(nrow(new_fit@data@y), 10)

  # Check parboot (SSE doesn't work)
  pb <- parboot(new_fit, statistic=function(x) x@AIC, nsim=2)
  expect_equal(pb@t.star[1,1], 109.8942, tol=1e-4)

  # Nonparboot
  npb <- nonparboot(new_fit, B=2)
  expect_equal(length(npb@bootstrapSamples), 2)
  expect_equal(npb@bootstrapSamples[[1]]@AIC, 136.29, tol=1e-4)
  expect_equal(numSites(npb@bootstrapSamples[[1]]@data), numSites(npb@data))
  v <- vcov(npb, method='nonparboot')
  expect_equal(nrow(v), length(coef(npb)))

  # getP
  gp <- getP(m1)
  expect_equal(dim(gp), dim(m1@data@y))
  expect_equal(as.vector(gp[5:6, 5:6]), c(0.4513,0.4513,0.7232,0.7232), tol=1e-4)

  # ranef
  expect_error(r <- ranef(m1))

  # Check error when random effect in formula
  expect_error(occuFP(~(1|dummy), ~1, ~1, data=umf1))
})
