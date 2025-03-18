context("occuComm fitting function")
skip_on_cran()

nsite <- 300
nocc <- 5
nsp <- 30

set.seed(123)
# Create a site by species covariate
x <- matrix(rnorm(nsite*nsp), nsite, nsp)

mu_0 <- 0
sd_0 <- 0.4
beta0 <- rnorm(nsp, mu_0, sd_0)

mu_x <- 1
sd_x <- 0.3
beta_x <- rnorm(nsp, mu_x, sd_x)

mu_a <- 0
sd_a <- 0.2
alpha0 <- rnorm(nsp, mu_a, sd_a)

ylist <- list()
z <- matrix(NA, nsite, nsp)
for (s in 1:nsp){
  psi <- plogis(beta0[s] + beta_x[s] * x[,s])
  z[,s] <- rbinom(nsite, 1, psi)

  p <- plogis(alpha0[s])

  y <- matrix(NA, nsite, nocc)
  for (m in 1:nsite){
    y[m,] <- rbinom(nocc, 1, p * z[m,s])
  }
  ylist[[s]] <- y
}
names(ylist) <- paste0("sp", sprintf("%02d", 1:nsp))
sc <- data.frame(a=factor(sample(letters[1:5], nsite, replace=TRUE)))
oc <- list(b = matrix(rnorm(nsite*nocc), nsite, nocc))

test_that("unmarkedFrameOccuComm can be constructed",{

  # Standard umf
  umf <- unmarkedFrameOccuComm(ylist, siteCovs=sc, obsCovs=oc)
  expect_is(umf, "unmarkedFrameOccuComm")
  expect_equal(numSites(umf), nsite)
  expect_equal(names(siteCovs(umf)), "a")
  expect_equal(names(obsCovs(umf)), "b")
  expect_null(umf@speciesCovs)
  expect_is(umf@ylist, "list")
  expect_equal(length(umf@ylist), nsp)

  # Provide y as array, gets converted to list
  yarray <- array(NA, c(nsite, nocc, nsp))
  for (i in 1:nsp){
    yarray[,,i] <- ylist[[i]]
  }
  umf2 <- unmarkedFrameOccuComm(yarray)
  expect_identical(umf2@ylist, umf@ylist)

  # Species covs handling
  # unmarked will reject all covariates with the wrong dimensions, e.g.
  spc_bad <- list(c = matrix(rnorm(nsite*nsp), nsite, nsp), 
                d = array(rnorm(nsite*nsp*nocc), c(nsite, nocc, nsp)),
                e = matrix(rnorm(nsite*nocc), nsite, nocc))
  expect_error(unmarkedFrameOccuComm(ylist, sc, speciesCovs=spc_bad),
               "Species covariate")

  # The one we created above works because it's M x S
  spc <- list(x = x,c = matrix(rnorm(nsite*nsp), nsite, nsp), 
                d = array(rnorm(nsite*nsp*nocc), c(nsite, nocc, nsp)),
                e = rnorm(nsp))
  umf3 <- unmarkedFrameOccuComm(ylist, sc, speciesCovs = spc)
  expect_identical(umf3@speciesCovs, spc)

  # Brackets work
  umf4 <- umf3[1:10,]
  expect_equal(numSites(umf4), 10)
  expect_equal(nrow(umf4@speciesCovs$x), 10)
  expect_equal(nrow(umf4@speciesCovs$c), 10)
  expect_equal(length(umf4@speciesCovs$e), nsp)

  pl <- plot(umf) # test plotting

  nul <- capture.output(summary(umf))
  expect_equal(nul[1], "unmarkedFrameOccuComm Object")
})

test_that("multispeciesFormula works", {
  
  f1 <- multispeciesFormula(~1~1, NULL)
  expect_equal(f1$formula,
               ~1 + (1 | species) ~ 1 + (1 | species))

  # Covariates
  f2 <- multispeciesFormula(~a~b, NULL)
  expect_equal(f2$formula,
               ~a + (1 + a || species) ~ b + (1 + b || species))

  # Additional random effect
  f3 <- multispeciesFormula(~a~b + (1|x), NULL)
  expect_equal(f3$formula,
               ~a + (1 + a || species) ~ b + (1 + b || species) + (1 | x_species))

  # Length S species covariate
  f4 <- multispeciesFormula(~a~b, list(b=rnorm(3)))
  expect_equal(f4$formula,
               ~a + (1 + a || species) ~ b + (1 | species))
})

test_that("occuComm can fit models and methods work", {
  
  spc <- list(x = x)
  spc$x[1,] <- NA

  # add some NAs
  ylist_na <- lapply(ylist, function(x){
                       x[2,1] <- NA
                       x[3,] <- NA
                       x
                })

  umf <- unmarkedFrameOccuComm(ylist_na, sc, speciesCovs = spc)
 
  # quick fit
  fit <- expect_no_warning(occuComm(~1~x, umf[1:20,]))
  expect_equivalent(coef(fit),
               c(0.1262, 1.3458, -0.00887), tol=1e-4)

  # Methods

  # summary
  nul <- capture.output(summary(fit))
  expect_equal(nul[2], "Call:")

  # predict
  pr <- predict(fit, 'state')
  expect_is(pr, "list")
  expect_equal(names(pr), names(ylist))
  expect_equal(dim(pr[[1]]), c(20, 4))
  expect_true(is.na(pr[[1]]$Predicted[1]))
  expect_equal(pr[[1]]$Predicted[3], 0.9423, tol=1e-4)
  
  pr <- predict(fit, 'det')
  expect_is(pr, "list")
  expect_equal(names(pr), names(ylist))
  expect_equal(dim(pr[[1]]), c(20*nocc, 4))
  expect_equal(pr[[1]]$Predicted[1], 0.5754, tol=1e-4)

  nd <- data.frame(x = c(0, 1))
  # no species
  expect_error(predict(fit, 'state', newdata=nd), "species")

  nd$species <- "sp01"
  pr <- predict(fit, 'state', newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equal(pr$Predicted[1], 0.6139, tol=1e-4)

  # sigma
  sig <- sigma(fit)
  expect_equal(nrow(sig), 3)

  # randomTerms
  rterm <- randomTerms(fit)
  expect_equal(nrow(rterm), 3 * length(ylist))
  expect_equal(rterm$Level, rep(names(ylist), 3))

  # Add the means
  rterm2 <- randomTerms(fit, addMean = TRUE)
  expect_equivalent(rterm2$Estimate, rterm$Estimate + rep(coef(fit), each=30))

  # richness
  rich <- richness(fit)
  expect_equal(length(rich), 20)
  rich_post <- richness(fit, posterior=TRUE)
  expect_equal(dim(rich_post@samples), c(20, 1, 100))

  # p
  gp <- getP(fit)
  expect_is(gp, "list")
  expect_equal(names(gp), names(ylist))
  expect_equal(dim(gp[[1]]), c(20, nocc))
  expect_equal(gp[[1]][1,1], 0.57535, tol=1e-4)

  # residuals
  res <- residuals(fit)
  expect_is(res, "list")
  expect_equal(names(res), names(ylist))
  expect_equal(dim(res[[1]]), c(20, nocc))
  expect_true(is.na(res[[1]][1,1]))
  expect_true(all(is.na(res[[1]][3,])))
  expect_equal(res[[1]][2,2], -0.3049, tol=1e-4)

  pdf(NULL)
  pl <- plot(fit) # test residual plot
  dev.off()

  # sse
  sse <- SSE(fit)
  expect_equivalent(sse, 460.1721, tol=1e-4)

  # update
  upd <- expect_no_warning(update(fit, formula=~1~1))
  expect_equal(length(coef(upd)), 2)
  upd <- expect_no_warning(update(fit, data=umf[1:10,]))
  expect_equal(numSites(upd@data), 10)

  # plotEffects
  # FIXME: doesn't work with species covs
  expect_error(plotEffects(fit, 'state', covariate='x'))

  # simulate
  s <- expect_warning(simulate(fit, nsim=2))
  expect_equal(length(s), 2)
  expect_equal(length(s[[1]]), nsp)
  expect_equal(dim(s[[1]][[1]]), c(20, nocc))
  expect_true(all(is.na(s[[1]][[1]][1,])))

  # parboot
  pb <- expect_warning(parboot(fit, nsim=2))

  # nonparboot
  npb <- expect_no_warning(nonparboot(fit, B=2))
  expect_equal(length(npb@bootstrapSamples), 2)
})

test_that("Additional random effect can be fitted", {
  
  set.seed(123)
  spc <- list(x=x)
  sc$new <- factor(sample(letters[1:10], nsite, replace=TRUE))
  umf <- unmarkedFrameOccuComm(ylist, sc, speciesCovs = spc)
 
  # quick fit
  fit <- occuComm(~1~x + (1|new), umf[1:20,])

  expect_equal(nrow(sigma(fit)), 4)
  expect_equivalent(length(coef(fit)), 3)

  # only random intercepts allowed
  expect_error(occuComm(~1~x + (x|new), umf[1:20,]),
               "random intercepts")

})

test_that("length S species covs can be fit", {
 
  set.seed(123)
  spc <- list(x = rnorm(nsp))
  umf <- unmarkedFrameOccuComm(ylist, sc, speciesCovs = spc)
 
  # quick fit
  fit <- occuComm(~1~x, umf[1:20,])

  # only 2 random effects, since x is length S
  expect_equal(expect_warning(nrow(sigma(fit))), 2)
  expect_equal(length(coef(fit)), 3) # still 3 coefficients

})
