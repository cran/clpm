\name{vcov.clpm}
\alias{vcov.clpm}
\title{
Calculate Variance-Covariance Matrix for a Fitted Constrained Linear Probablity Model Object
}
\description{
Returns the variance-covariance matrix of the main parameters of a fitted constrained linear probability model object. 
}

\usage{
\method{vcov}{clpm}(object, \ldots)
}

\arguments{
  \item{object}{a fitted model object, typically.}
  \item{\ldots}{additional arguments for method functions.}
}

\details{
See \link{vcov.lm}.
}

\value{
The function returns the same object returned by \link{vcov.lm}.
}

\author{Andrea Beci \email{andreabeci08@gmail.com}, Paolo Frumento \email{paolo.frumento@unipi.it}}
\keyword{methods}