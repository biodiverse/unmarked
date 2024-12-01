\name{unmarkedFrameOccuComm}
\alias{unmarkedFrameOccuComm}

\title{Organize data for the community occupancy model fit by occuComm}

\usage{unmarkedFrameOccuComm(y, siteCovs = NULL, obsCovs = NULL, speciesCovs = NULL)}

\description{Organize detection/non-detection data and covariates for use with
  the community occupancy model}

\arguments{
  \item{y}{The detection-nondetection data for multiple species. This data may
   be in one of two forms: either a named list of S M x J matrices (M = sites, 
   J = occasions, S = species) or an M x J x S array.
  }
  \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
    site level. This should have M rows and one column per covariate.}
  \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
    covariates that vary within sites, or a \code{\link{data.frame}} with 
    M x J rows in site-major order.}
  \item{speciesCovs}{Covariates that also vary by species. This must be provided
    as a named list. Each list element can have one of three possible dimensions:
    a vector of length S (e.g., mean species body mass); a matrix M x S
    (covariates that vary by site and species); or an array M x J x S
    (covariates that vary by site, occasion, and species).
  }
}

\value{an object of class unmarkedFrameOccuComm}

\seealso{\code{\link{occuComm}}}

\examples{

# Simulate some multispecies data
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
names(ylist) <- paste0("sp", sprintf("\%02d", 1:nsp))
sc <- data.frame(a=factor(sample(letters[1:5], nsite, replace=TRUE)))

# Species covs need to be a list of named elements
spc <- list(x = x)

umf <- unmarkedFrameOccuComm(ylist, sc, speciesCovs = spc)
summary(umf)

}
