\name{openmp}
\alias{openmp}
\title{Control number of OpenMP threads used by TMB engine}
\description{
This function sets the number of parallel threads that will be used by OpenMP if it
is available on your system. Only applies to models fit with the TMB engine.
}

\usage{
openmp(threads = 1)
}

\arguments{
  \item{threads}{The number of threads to use}
}

\value{The number of threads.}
