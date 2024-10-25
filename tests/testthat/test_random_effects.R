context("random effects tools")
skip_on_cran()

set.seed(123)
M <- 10
J <- 3

sc <- data.frame(x1 = rnorm(M),
                 x2 = factor(sample(letters[1:2], M, replace=TRUE)),
                 x3 = factor(sample(letters[3:4], M, replace=TRUE)),
                 g = factor(sample(letters[5:7], M, replace=TRUE)))

test_that("random factor slopes work", {

  test1 <- get_Z(~x1 + x2 + (x1 + x2 || g), sc)
  expect_is(test1, "dgCMatrix")

  test1 <- as.matrix(test1)
  expect_equal(colnames(test1), rep(letters[5:7], 4))
  
  expect_equal(test1, structure(c(0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1,
0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1.55870831414912,
0, 0.129287735160946, 1.71506498688328, 0, 0, 0, -0.445661970099958,
0, -0.23017748948328, 0, 0, 0, 0, 0.460916205989202, 0, 0, 0,
-0.560475646552213, 0, 0, 0.070508391424576, 0, 0, 0, -1.26506123460653,
-0.686852851893526, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 1, 0, 0), dim = c(10L, 12L), dimnames = list(c("1", "2",
"3", "4", "5", "6", "7", "8", "9", "10"), c("e", "f", "g", "e",
"f", "g", "e", "f", "g", "e", "f", "g"))))

  sig <- sigma_names(~x1 + x2 + (x1 + x2 || g), sc)
  expect_equal(sig,
               c("1|g", "x1|g", "x2a|g", "x2b|g"))

  nr <- get_nrandom(~x1 + x2 + (x1 + x2 || g), sc)
  expect_equal(length(nr), length(sig))
  expect_equal(sum(nr), ncol(test1))

  test2 <- as.matrix(get_Z(~0 + x2:x3 + (0 + x2:x3 | g), sc))
  sig <- sigma_names(~0 + x2:x3 + (0 + x2:x3 | g), sc)
  expect_equal(sig, c("x2a:x3c|g", "x2b:x3c|g", "x2a:x3d|g", "x2b:x3d|g"))

  expect_equal(test2, structure(c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), dim = c(10L,
12L), dimnames = list(c("1", "2", "3", "4", "5", "6", "7", "8",
"9", "10"), c("e", "f", "g", "e", "f", "g", "e", "f", "g", "e",
"f", "g"))))
  
  nr <- get_nrandom(~0 + x2:x3 + (0 + x2:x3 | g), sc)
  expect_equal(length(nr), length(sig))
  expect_equal(sum(nr), ncol(test2))
})

test_that("random factor slope estimates work", {
skip_on_ci()

M <- 10000
J <- 3

set.seed(123)
sc <- data.frame(x1 = rnorm(M),
                 x2 = factor(sample(letters[1:2], M, replace=TRUE)),
                 x3 = factor(sample(letters[3:4], M, replace=TRUE)),
                 g = factor(sample(letters[5:26], M, replace=TRUE)))

X <- model.matrix(~x1 + x2, sc)

b <- c(-0.3, 0.3, 0.5)

b2 <- c(-0.3, 0.3, 0, 0.5)

re <- lapply(1:length(b2), function(i){
  rnorm(length(levels(sc$g)), 0, 0.2)
})

Z <- as.matrix(unmarked:::get_Z(~x1 + x2 + (x1 + x2 ||g), sc))

Xall <- cbind(X, Z)

ball <- c(b, unlist(re))

psi <- plogis(Xall %*% ball)
z <- rbinom(M, 1, psi)

y <- matrix(NA, M, J)
for (i in 1:M){
  y[i,] <- rbinom(3, 1, 0.4) * z[i]
}

umf <- unmarkedFrameOccu(y=y, siteCovs=sc)

fit <- occu(~1~x1 + x2 + (x1 + x2 || g), umf)

expect_equivalent(coef(fit),
             c(b, qlogis(0.4)), tol=0.1)

expect_equivalent(sigma(fit)$sigma,
                  c(0.2865,0.1859,0.3004,0.1928), tol=1e-4)

re_est <- randomTerms(fit)
expect_equal(re_est$Level, rep(levels(sc$g), 4))

pr <- predict(fit, 'state', level=NULL)$Predicted

expect_true(cor(psi, pr) > 0.95)

# Two grouping variables
fit2 <- occu(~1~x1 + x2 + (x1 + x2 || g) + (1|x3), umf[1:100,])

expect_equal(nrow(sigma(fit2)), 5)

re <- randomTerms(fit2)
expect_equal(re$Groups[89:90], c("x3", "x3"))
})
