\name{clpm-package}
\alias{clpm-package}
\docType{package}
\title{
Constrained Estimation of the Linear Probability Model
}
\description{
Implementation of a method to estimate the linear probability model, ensuring that the predicted probabilities are inside the (0,1) interval. The same method can be applied to any variable for which predictions between 0 and 1 are required.
}
\details{
\tabular{ll}{
Package: \tab clpm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2023-04-30\cr
License: \tab GPL-2\cr
}
The function \code{\link{clpm}}
}
\author{
Andrea Beci, Paolo Frumento

Maintainer: Andrea Beci <andreabeci08@gmail.com>
}
\keyword{ package }
\examples{
x <- runif(100)
y <- rbinom(100, 1, x)
clpm(y~x)
}




