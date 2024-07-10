\name{unmarkedPower-methods}
\alias{unmarkedPower-methods}
\alias{unmarkedPower-class}
\alias{unmarkedPower_old-class}
\alias{show,unmarkedPower-method}
\alias{summary,unmarkedPower-method}
\alias{update,unmarkedPower-method}
\alias{show,unmarkedPower_old-method}
\alias{summary,unmarkedPower_old-method}
\alias{update,unmarkedPower_old-method}
\alias{plot,unmarkedPower,missing-method}

\title{Methods for unmarkedPower objects}

\description{Various functions to summarize unmarkedPower objects}

\usage{
\S4method{show}{unmarkedPower}(object)
\S4method{summary}{unmarkedPower}(object, alpha, showIntercepts = FALSE, ...)
\S4method{plot}{unmarkedPower,missing}(x, y, alpha, showIntercepts = FALSE, ...)
}

\arguments{
  \item{object,x}{An object of class \code{unmarkedPower} created with the 
    \code{powerAnalysis} function}
  \item{alpha}{Desired Type I error rate. If not provided, defaults to the value
    specified when calling \code{powerAnalysis}.}
  \item{showIntercepts}{Show intercepts output? This is 
    rarely useful.}
  \item{y}{Not currently used.}
  \item{...}{Not currently used.}
}

\value{
  For \code{show} and \code{summary}, summary output is printed to the console.
  For \code{plot}, a visualization of the summary output is created.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{powerAnalysis}}
}
