\name{clpm.control}
\alias{clpm.control}
\title{
Control Parameters for Constrained Estimation of Linear Probability Model
}
\description{
This function generates default values for the tolerance (\code{tol}) and maximum number of iterations (\code{maxit}) for \link{clpm}, as well as an option to print the progress.
}

\usage{
clpm.control(tol = 1e-10, maxit, trace = FALSE)
}

\arguments{
  \item{tol}{tolerance value for the steps of the algorithm.}
  \item{maxit}{maximum number of attempts of the algorithm.}
  \item{trace}{logical: should the progress be printed on screen?}
}

\value{
The function performs a sanity check and returns its arguments.
}

\author{Andrea Beci \email{andreabeci08@gmail.com}, Paolo Frumento \email{paolo.frumento@unipi.it}}
\keyword{methods}

