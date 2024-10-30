\name{richness}
\alias{richness}
\alias{richness,unmarkedFitOccuComm-method}

\title{Estimate posterior distributions of site richness}
\description{
 Estimate posterior distributions of site richness from a fitted community occupancy
 model using empirical Bayes methods. This function uses \code{\link{ranef}}
 and \code{\link{posteriorSamples}} internally.
}

\usage{
\S4method{richness}{unmarkedFitOccuComm}(object, nsims=100, posterior=FALSE)
}

\arguments{
  \item{object}{A \code{unmarkedFitOccuComm} object.}
  \item{nsims}{Number of posterior samples to generate.} 
  \item{posterior}{Logical. If \code{FALSE} mean estimates for each site are
    returned; if \code{TRUE} the entire posterior for each site is returned.}
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{occuComm}}, \code{\link{ranef}}, \code{\link{posteriorSamples}}}
