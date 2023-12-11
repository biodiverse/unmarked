\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,unmarkedFit-method}
\alias{predict,unmarkedFitOccuFP-method}
\alias{predict,unmarkedFitOccuMulti-method}
\alias{predict,unmarkedFitOccuMS-method}
\alias{predict,unmarkedFitOccuTTD-method}
\alias{predict,unmarkedFitNmixTTD-method}
\alias{predict,unmarkedFitPCount-method}
\alias{predict,unmarkedFitColExt-method}
\alias{predict,unmarkedFitGMM-method}
\alias{predict,unmarkedFitGDS-method}
\alias{predict,unmarkedFitPCO-method}
\alias{predict,unmarkedFitDSO-method}
\alias{predict,unmarkedFitGDR-method}
\alias{predict,unmarkedFitList-method}
\alias{predict,unmarkedRanef-method}
\title{ Methods for Function predict in Package `unmarked' }

\usage{
predict(object, type, newdata, 
        backTransform = TRUE, na.rm = TRUE,
        appendData = FALSE, level = 0.95, 
        re.form = NULL, ...)
}

\description{
These methods return predicted values from fitted model objects.
}

\arguments{

  \item{object}{An \linkS4class{unmarkedFit} object or an \linkS4class{unmarkedFitList} object. See methods.}
  
  \item{type}{A character string that depends upon the fitted models. You can see the possible types with \code{names(umfit)}, with \code{umfit} your \code{unmarkedFit} object.}
  
  \item{newdata}{Facultative: if no newdata, uses the actual data (\code{\link{getData}(umfit)}). When specified, \code{newdata} must be \code{\link{unmarkedFrame}}, \code{\link{data.frame}}, \code{\link{RasterLayer}}, \code{\link{RasterStack}}, or \code{\link{SpatRaster}}. See example.}
  
  \item{backTransform}{Logical specifying whether to back-transform parameters (\code{backTransform=TRUE}) or not (\code{backTransform=FALSE}).}

  \item{na.rm}{Logical specifying whether to remove NAs (\code{na.rm=TRUE}) or not (\code{na.rm=FALSE}).}

  \item{appendData}{Logical specifying whether to append the data (covariates) used for the prediction (\code{appendData=TRUE}) or not (\code{appendData=FALSE}).}

  \item{level}{Numeric specifying the confidence interval level.}
  
  \item{re.form}{}
  
  \item{\dots}{Unused.}
  
}


\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFit")}}{
"type" depends upon the fitted models. You can see the possible types with \code{names(umfit), with \code{umfit} your \code{unmarkedFit} object.}
}
\item{\code{signature(object = "unmarkedFitList")}}{
"type" depends upon the fitted models. You can see the possible types with \code{names(fl), with \code{fl} your \code{unmarkedFitList} object.}
}
\item{\code{signature(object = "unmarkedFitOccu")}}{
"type" must be either 'state' or 'det'.
}
\item{\code{signature(object = "unmarkedFitColExt")}}{
"type" must be 'psi', 'col', 'ext', or 'det'.
}
\item{\code{signature(object = "unmarkedFitGMM")}}{
"type" must be 'lambda', 'psi', 'det'
}
\item{\code{signature(object = "unmarkedRanef")}}{
Use this method to generate the empirical Bayes posterior predictive distribution
for functions of the random variables (latent abundance or occurrence).

In addition to the output object from \code{ranef}, you must also supply a 
custom function to argument \code{func}. The function must take as input a matrix
with dimensions M x T, where M is the number of sites and T is the number of 
primary periods (T=1 for single-season models). The output of this function should 
be a vector or matrix containing the derived parameters of interest. 

You may also manually set the number of draws from the posterior predictive
distribution with argument \code{nsims}; the default is 100.

The output of \code{predict} will be a vector or array with one more dimension
than the output of the function supplied \code{func}, corresponding to the number
of draws requested \code{nsims}. For example, if \code{func}
outputs a scalar, the output of \code{predict} will be a vector with length
equal to \code{nsims}. If \code{func} outputs a 3x2 matrix, the output of
\code{predict} will be an array with dimensions 3x2x\code{nsims}.
See \code{\link{ranef}} for an example.

Alternatively, you can use the \code{\link{posteriorSamples}} function on the
\code{ranef} output object to obtain the full posterior predictive distribution. 
This is useful if you are having trouble designing your custom function or if 
you want to obtain multiple different derived parameters from the same posterior
predictive distribution.
}
\item{\code{signature(object = "unmarkedFitCOP")}}{
"type" must be either 'psi' or 'lambda'.
}

}}

\examples{
set.seed(123)
options(max.print = 100)

# Simulate an unmarkedFrame with count data in 100 sites for 3 sampling occasions
(umf <- simulate(
  "COP",
  formulas = list(psi =  ~ habitat, lambda =  ~ rain),
  coefs = list(
    psi = c(intercept = qlogis(.2), habitatForest = .5, habitatGrassland = .8),
    lambda = c(intercept = log(1), rain = -1)
  ),
  design = list(M = 100, J = 3),
  guide = list(habitat = factor(levels = c("City", "Forest", "Grassland")))
))

# Fit the model
(umfit <- occuCOP(umf, psiformula = ~ habitat, lambdaformula =  ~ rain, L1 = TRUE))

# See the possible types for this unmarkedFit object
names(umfit)

# Back-transform the occupancy probability
# = predict with the data furnished to the model
predict(object = umfit, type = "psi", appendData = TRUE)

# Because occupancy only depends on the habitat, 
# we have the same information in a more concise form by specifying newdata
# as the 3 habitat types in our covariates data.
predict(
  object = umfit,
  type = "psi",
  newdata = data.frame("habitat" = c("Forest", "Grassland", "City")),
  appendData = TRUE
)

# We can look at the 90% confidence interval (instead of the default 95% CI)
predict(
  object = umfit,
  type = "psi",
  newdata = data.frame("habitat" = c("Forest", "Grassland", "City")),
  level = 0.90,
  appendData = TRUE
)

# Back-transform the detection rate
predict(object = umfit, type = "lambda", appendData = TRUE)

## This is not easily readable. We can show the results in a clearer way, by:
##  - adding the site and observation
##  - printing only the wind covariate used to get the predicted lambda
cbind(
  data.frame(
    "site" = rep(1:numSites(umf), each = obsNum(umf)),
    "obs" = rep(1:obsNum(umf), times = numSites(umf)),
    "rain" = getData(umfit)@obsCovs
  ),
  predict(umfit, "lambda", appendData = FALSE)
)

# We can predict the detection rate with new detection covariates
predict(umfit, "lambda", 
        newdata = data.frame("rain" = c(0, 1, 2)),
        appendData = TRUE)
}

\keyword{methods}
