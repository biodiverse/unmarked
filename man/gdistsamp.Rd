\name{gdistsamp}
\alias{gdistsamp}
\title{
  Fit the generalized distance sampling model of Chandler et al. (2011).
}
\description{
  Extends the distance sampling model of Royle et al. (2004) to estimate
  the probability of being available for detection. Also allows
  abundance to be modeled using the negative binomial and zero-inflated
  Poisson distributions. 
}
\usage{
gdistsamp(lambdaformula, phiformula, pformula, data, keyfun =
c("halfnorm", "exp", "hazard", "uniform"), output = c("abund",
"density"), unitsOut = c("ha", "kmsq"), mixture = c("P", "NB", "ZIP"), K,
starts, method = "BFGS", se = TRUE, engine=c("C","R"), rel.tol=1e-4, threads=1, ...)
}
\arguments{
  \item{lambdaformula}{
      A right-hand side formula describing the abundance covariates.
}
  \item{phiformula}{
      A right-hand side formula describing the availability covariates.
}
  \item{pformula}{
      A right-hand side formula describing the detection function covariates.
}
  \item{data}{
      An object of class \code{unmarkedFrameGDS}
}
  \item{keyfun}{
      One of the following detection functions: "halfnorm", "hazard", "exp",
          or "uniform." See details.
}
  \item{output}{
      Model either "density" or "abund"
}
  \item{unitsOut}{
      Units of density. Either "ha" or "kmsq" for hectares and square
          kilometers, respectively.
}
  \item{mixture}{
      Either "P", "NB", or "ZIP" for the Poisson, negative binomial, or
      zero-inflated Poisson models of abundance.
}
  \item{K}{
      An integer value specifying the upper bound used in the integration.
}
  \item{starts}{
      A numeric vector of starting values for the model parameters.
}
  \item{method}{
      Optimization method used by \code{\link{optim}}.
}
  \item{se}{
      logical specifying whether or not to compute standard errors.
}
  \item{engine}{
      Either "C" to use fast C++ code or "R" to use native R code during the 
      optimization.
}
\item{rel.tol}{relative accuracy for the integration of the detection function.
    See \link{integrate}. You might try adjusting this if you get an error
    message related to the integral. Alternatively, try providing
    different starting values.}

\item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled.} 

\item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}

}

\details{

  Extends the model of Royle et al. (2004) by estimating the probability
  of being available for detection \eqn{\phi}{phi}. To estimate this
  additional parameter, replicate distance sampling data must be
  collected at each transect. Thus the data are collected at i = 1, 2,
  ..., R transects on t = 1, 2, ..., T occassions. As with the model of
  Royle et al. (2004), the detections must be binned into distance
  classes. These data must be formatted in a matrix with R rows, and JT
  columns where J is the number of distance classses. See
  \code{\link{unmarkedFrameGDS}} for more information about data
  formatting.

  The definition of availability depends on the context. The model is
  \deqn{M_i \sim \text{Pois}(\lambda)}{M(i)~Pois(lambda)}
  \deqn{N_{i,t} \sim \text{Bin}(M_i, \phi)}{N(i,t)~Bin(M(i), phi)}
  \deqn{y_{i,1,t}, \dots, y_{i,J,t} \sim \text{Multinomial}(N_{i,t},
  \pi_{i,1,t}, \dots, \pi_{i,J,t})}{y(i,1,t), ..., y(i,J,t) ~
  Multinomial(N(i,t), pi(i,1,t), ..., pi(i,J,t))} 

  If there is no movement, then \eqn{M_i}{M(i)} is local abundance, and
  \eqn{N_{i,t}}{N(i,t)} is the number of individuals that are available
  to be detected. In this case, \eqn{\phi=g_0}{phi=g(0)}. Animals might
  be missed on the transect line because they are difficult to see or
  detected. This relaxes the assumption of conventional distance
  sampling that \eqn{g_0=1}{g(0)=1}. 

  However, when there is movement in the form of temporary emigration,
  local abundance is \eqn{N_{i,t}}{N(i,t)};  it's the fraction of
  \eqn{M_i}{M(i)} that are on the plot at time t. In this case,
  \eqn{\phi}{phi} is the temporary emigration parameter, and we need to
  assume that \eqn{g_0=1}{g(0)=1} in order to interpret
  \eqn{N_{i,t}}{N(i,t)} as local abundance. See Chandler et al. (2011)
  for an analysis of the model under this form of temporary emigration.

  If there is movement and \eqn{g_0<1}{g(0)<1} then it
  isn't possible to estimate local abundance at time t. In this case,
  \eqn{M_i}{M(i)} would be the total number of individuals that ever use
  plot i (the super-population), and \eqn{N_{i,t}}{N(i,t)} would be the
  number available to be detected at time t. Since a fraction of the
  unavailable individuals could be off the plot, and another fraction
  could be on the plot, it isn't possible to infer local abundance and
  density during occasion t.   

}

\note{
  If you aren't interested in estimating \eqn{\phi}{phi}, but you want
  to use the negative binomial or ZIP distributions, set numPrimary=1
  when formatting the data. 
}

\value{
  An object of class unmarkedFitGDS.
}

\references{

  Royle, J. A., D. K. Dawson, and S. Bates. 2004. Modeling abundance
    effects in distance sampling. \emph{Ecology} 85:1591-1597.

  Chandler, R. B, J. A. Royle, and D. I. King. 2011. Inference about
    density and temporary emigration in unmarked
    populations. \emph{Ecology}  92:1429--1435. 
}

\author{
  Richard Chandler \email{rbchan@uga.edu}
}

\note{
  You cannot use obsCovs, but you can use yearlySiteCovs (a confusing
  name since this model isn't for multi-year data. It's just a hold-over
  from the colext methods of formatting data upon which it is based.) 
}

\seealso{
  \code{\link{distsamp}}
}

\examples{


# Simulate some line-transect data

set.seed(36837)

R <- 50 # number of transects
T <- 5  # number of replicates
strip.width <- 50
transect.length <- 100
breaks <- seq(0, 50, by=10)

lambda <- 5 # Abundance
phi <- 0.6  # Availability
sigma <- 30 # Half-normal shape parameter

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda) # Individuals within the 1-ha strip
    for(t in 1:T) {
        # Distances from point
        d <- runif(M, 0, strip.width)
        # Detection process
        if(length(d)) {
            cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
            d <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R) # convert array to matrix

# Organize data
umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
    dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)


# Fit the model
m1 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=50)

summary(m1)

predict(m1, type = 'lambda')[1,]
predict(m1, type = 'phi')[1,]
predict(m1, type = 'det')[1,]

\dontrun{
# Empirical Bayes estimates of abundance at each site
re <- ranef(m1)
plot(re, layout=c(10,5), xlim=c(-1, 20))
}

}
\keyword{ models }

