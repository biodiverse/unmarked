\name{unmarkedPowerList}
\alias{unmarkedPowerList}
\alias{unmarkedPowerList,unmarkedPower-method}
\alias{unmarkedPowerList-class}
\alias{unmarkedPowerList-methods}
\alias{show,unmarkedPowerList-method}
\alias{summary,unmarkedPowerList-method}
\alias{plot,unmarkedPowerList,ANY-method}

\alias{unmarkedPowerList_old-class}
\alias{unmarkedPowerList,list-method}
\alias{unmarkedPowerList,unmarkedFit-method}
\alias{show,unmarkedPowerList_old-method}
\alias{summary,unmarkedPowerList_old-method}
\alias{plot,unmarkedPowerList_old,ANY-method}

\title{Summarize a series of unmarked power analyses}

\description{
  A list of power analyses created with \code{powerAnalysis} can be combined 
  using \code{unmarkedPowerList}, allowing comparison e.g. between different
  study designs/sample sizes. A series of methods for \code{unmarkedPowerList}
  objects are available including a \code{plot} method.
}

\usage{
\S4method{unmarkedPowerList}{unmarkedPower}(object, ...)
\S4method{show}{unmarkedPowerList}(object)
\S4method{summary}{unmarkedPowerList}(object, showIntercepts = FALSE, ...)
\S4method{plot}{unmarkedPowerList,ANY}(x, power=NULL, param=NULL, ...)
}

\arguments{
  \item{object,x}{For \code{unmarkedPowerList}, an \code{unmarkedPower} object.
    For \code{show}, \code{summary}, \code{plot}, an \code{unmarkedPowerList} object.
  }
  \item{showIntercepts}{Show intercepts output? This is rarely useful.}
  \item{power}{If specified, adds a dotted line to the plot at this target power value.}
  \item{param}{When plotting, the model parameter to plot power vs. sample size for.
    By default this is the first parameter.}
  \item{...}{For \code{unmarkedPowerList}, other \code{unmarkedPower} objects
    to combine into the list.}
} 

\value{A \code{unmarkedPowerList} object, a summary of the object in the console,
      or a summary plot, depending on the method}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{powerAnalysis}}
}

\examples{

\dontrun{

# Build unmarkedFrame
umf <- unmarkedFrameOccu(y = matrix(NA, 300, 8),
                         siteCovs = data.frame(elev=rnorm(300)))

# Run power analyses
cf <- list(state = c(0, -0.4), det = 0)
pa1 <- powerAnalysis(umf, model=occu, formula=~1~elev, effects=cf)
pa2 <- powerAnalysis(umf[1:100,], model=occu, formula=~1~elev, effects=cf)

# Combine them into a list
(pl <- unmarkedPowerList(pa1, pa2))

# Look at summary plot for elev effect
plot(pl, power=0.8, param='elev')

}
}
