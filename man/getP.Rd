\name{getP}
\alias{getP}
\alias{getP,unmarkedFit-method}

\title{Get detection probability matrix from a fitted model}

\description{
  Calculates detection probability for each site and occasion and returns the values
  as a matrix with (generally) the same dimensions as the y-matrix.
}

\usage{
\S4method{getP}{unmarkedFit}(object, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{...}{Other arguments, currently ignored.}
}

\value{A matrix of estimated detection probabilities with the 
  same dimensions as the \code{y} matrix (M x J or M x JT). For some models with 
  multiple observation processes or multiple species, such as \code{gdistremoval}
  or \code{occuMulti}, the output may be a named list of matrices.
}
