\name{predict.clpm}
\alias{predict.clpm}
\title{
Predict method for Constrained Linear Probability Model
}
\description{
Predicted values based on constrained linear probability model object.
}

\usage{
\method{predict}{clpm}(object, newdata, se.fit = FALSE, \ldots)
}

\arguments{
  \item{object}{Object of class inheriting from \code{"clpm"}}
  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{se.fit}{A switch indicating if standard errors are required.}
  \item{\ldots}{further arguments passed to or from other methods.}
}

\details{
See \link{predict.lm}.
}

\value{
The function returns the same object returned by \link{predict.lm}.
}

\author{Andrea Beci \email{andreabeci08@gmail.com}, Paolo Frumento \email{paolo.frumento@unipi.it}}
\keyword{methods}