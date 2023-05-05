\name{clpm}
\alias{clpm}
\title{
Constrained Estimation of the Linear Probability Model
}
\description{
\code{clpm} is used to fit the linear probability model while ensuring that the predicted probabilities are in the (0,1) interval. The function can also be applied to any variable for which predictions between 0 and 1 are required.
}
\usage{
clpm(formula , data, subset, na.action, weights, contrasts = NULL,
lambda = NULL, control = clpm.control(), \ldots)
}
\arguments{
  \item{formula}{a two-sided formula of the form \code{y ~ x1 + x2 + \ldots}: a symbolic description of the linear probability model.
  The \code{y} argument must be a variable (binary or continous) for which we want predictions to be inside the (0,1) interval. The model specification is exactly as in \link{lm}.
  }
  \item{data}{
    an optional data frame, list or environment (or object coercible by \link{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{clpm} is called.
  }
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs. See \link{lm} for details.}
  \item{weights}{an optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector. See \link{lm}.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \link{model.matrix.default}.}
  \item{lambda}{a tuning constant that defines how important it is to obtain predictions in the (0,1) interval. If \code{lambda} is too small, the constraints may not be respected. On the other hand, if \code{lambda} is too large, the objective function might lose its convexity. If no value is supplied,
  an optimal value will be selected iteratively.}
  \item{control}{see \link{clpm.control}.}
  \item{\ldots}{for future arguments.}
}
\details{For more details, see \link{lm}.}
\value{
\code{clpm} returns an object of class "\code{clpm}".

The functions \code{summary} and \code{predict} are used to obtain and print a summary and estimate model predictions.

An object of class \dQuote{\code{clpm}}, a list containing the following items:
\item{coefficients}{a named vector of coefficients.}
\item{covar}{the estimated variance-covariance matrix.}
\item{residuals}{the residuals, that is the response minus the fitted values.}
\item{rank}{the numeric rank of the fitted linear model.}
\item{fitted.values}{the fitted values, that represent conditional means or, for a binary response, conditional probabilities.}
\item{weights}{	(only for weighted fits) the specified weights.}
\item{df.residuals}{the residual degrees of freedom.}
\item{obj.function}{the value of the minimized loss function.}
\item{gradient}{the value of the gradient.}
\item{convergence}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{control}{the values from \link{clpm.control}.}
\item{lambda}{the lambda value applied for model estimation.}
\item{contrasts}{(only where relevant) the contrasts used.}
\item{xlevels}{(only where relevant) a record of the levels of the factors used in fitting.}
\item{call}{the matched call.}
\item{terms}{the \link{terms} object used.}
\item{model}{if requested (the default), the model frame used.}
}
\author{
Andrea Beci \email{andreabeci08@gmail.com}, Paolo Frumento \email{paolo.frumento@unipi.it}
}
\examples{
x <- runif(100)
y <- rbinom(100, 1, x)
fit <- clpm(y~x)

summary(fit)
predict(fit)
}
\keyword{models}
\keyword{regression}
