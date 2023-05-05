\name{summary.clpm}
\alias{summary.clpm}
\title{
Summary After Constrained Estimation of Linear Probability Model
}
\description{
Summary of an object of class \code{'clpm'}.
}

\usage{
\method{summary}{clpm}(object, correlation = FALSE, symbolic.cor = FALSE, \ldots)
}

\arguments{
  \item{object}{an object of class \code{"lm"}, usually, a result of a call to \link{clpm}.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form (see \link{symnum}) rather than as numbers.}
  \item{\ldots}{further arguments passed to or from other methods.}
}

\details{
See \link{summary.lm}.
}

\value{
The function returns the same object returned by \link{summary.lm}.
}

\author{Andrea Beci \email{andreabeci08@gmail.com}, Paolo Frumento \email{paolo.frumento@unipi.it}}
\keyword{methods}