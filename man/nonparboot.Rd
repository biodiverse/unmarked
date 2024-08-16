\name{nonparboot}
\alias{nonparboot}
\alias{nonparboot,unmarkedFit-method}

\title{Get non-parametric bootstrap samples from an \code{unmarked} model}
\description{
  Resamples sites (with replacement) from an \code{unmarkedFit} object and
  re-fits the model, to generate new non-parametric bootstrap samples of estimates.
  These samples can then be used by \code{vcov} in order to get bootstrap-based
  estimates of standard errors.
}

\usage{
\S4method{nonparboot}{unmarkedFit}(object, B = 1, keepOldSamples = TRUE, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{B}{The number of bootstrap samples to get}
  \item{keepOldSamples}{Should previously generated bootstrap samples 
   be retained?}
  \item{...}{Other arguments, currently ignored} 
}

\details{
  Calling \code{nonparboot} on an unmarkedFit returns the original
  unmarkedFit, with the bootstrap samples added on.
  Then subsequent calls to \code{\link{vcov}} with the argument
  \code{method="nonparboot"} will use these bootstrap samples.
  Additionally, standard errors of derived estimates from either
  \code{\link{linearComb}} or \code{\link{backTransform}} can be
  instructed to use bootstrap samples by providing the argument
  \code{method = "nonparboot"}.

  For all fitting functions, only sites are re-sampled. In previous versions of 
  \code{unmarked} both sites and occasions
  of \code{\link{occu}} and \code{\link{occuRN}} fits were re-sampled, but now
  only sites are re-sampled for these fit types to match the other fit types.

  Note: by default, if the model fit to a new resampled dataset results in an
  error or a convergence failure, \code{nonparboot} will continually 
  re-sample and re-fit until this is no longer the case. This could
  result in the function hanging.
}

\examples{
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), type = "removal")
(fm <- multinomPois(~ 1 ~ ufc + trba, ovenFrame))
fm <- nonparboot(fm, B = 20) # should use larger B in real life.
vcov(fm, method = "hessian")
vcov(fm, method = "nonparboot")
avg.abundance <- backTransform(linearComb(fm, type = "state", coefficients = c(1, 0, 0)))

## Bootstrap sample information propagates through to derived quantities.
vcov(avg.abundance, method = "hessian")
vcov(avg.abundance, method = "nonparboot")
SE(avg.abundance, method = "nonparboot")
}
\keyword{methods}
