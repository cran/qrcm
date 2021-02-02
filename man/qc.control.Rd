\name{qc.control}
\alias{qc.control}
\title{
Estimate Non-Crossing Quantile Functions
}
\description{
This function generates a list of arguments to be used as operational parameters for \code{remove.qc} within a call to \code{\link{iqr}}. Additionally, this R documentation page contains a short description of the algorithm, which is presented in details in Sottile and Frumento (2021). 
}
\usage{
qc.control(maxTry = 25, trace = FALSE, lambda = NULL)
}
\arguments{
  \item{maxTry}{
  maximum number of attempts of the algorithm.
}
  \item{trace}{
  logical: should the progress be printed on screen?
}
  \item{lambda}{
  an optional positive scalar to be used as tuning parameter (see \dQuote{Details}). By default, \kbd{lambda = NULL}.
}
}
\details{
Quantile crossing occurs when the first derivative of the estimated quantile function is negative at some value of \eqn{p}. The argument \code{remove.qc} of the \command{iqr} function can be used to eliminate quantile crossing. 

The algorithm proceeds as follows. A penalization that reflects the severity of crossing is added to the loss function. The weight of the penalty term is determined by a tuning parameter \eqn{\lambda}. If \eqn{\lambda} is too small, the penalization has no effect. However, if \eqn{\lambda} is too large, the objective function may lose its convexity, causing a malfunctioning of the algorithm. In general, the value of \eqn{\lambda} is \emph{not} user-defined. The algorithm starts with an initial guess for the tuning parameter, and proceeds adaptively until it finds a suitable value. The maximum number of iterations is determined by the \kbd{maxTry} argument of this function (default \kbd{maxTry = 25}). The algorithm stops automatically when the \kbd{crossIndex} of the model (see \code{\link{diagnose.qc}}) is zero, or when no further progress is possible.

It is possible to supply a user-defined value of \eqn{\lambda}, e.g., \kbd{lambda = 7.5}. If this happens, the model is estimated \strong{once}, using the requested \kbd{lambda}, while the \kbd{maxTry} argument is ignored.


This method allows for censored or truncated data, that are supported by \command{iqr}. Full details are provided in Sottile and Frumento (2021).
}
\value{
The function performs a sanity check and returns its arguments.
}
\references{
Sottile, G., and Frumento, P. (2021). \emph{Parametric estimation of non-crossing quantile functions}. Statistical Modelling [forthcoming].
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\note{
Occasionally, the loss of the penalized model is smaller than that of the unconstrained fit. This is either an artifact due to numerical approximations or lack of convergence, or is explained by the fact that, if the quantile function is ill-defined, so is the loss function of the model. With censored or truncated data, however, it can also be explained by the fact that the \kbd{obj.function} of the model is \strong{NOT} the function being minimized (see note 3 in the documentation of \command{iqr}).
}


\seealso{
\code{\link{iqr}}, \code{\link{diagnose.qc}}.
}
\examples{
\donttest{
 # Using simulated data

 set.seed(1111)
 n <- 1000
 x1 <- runif(n,0,3)
 x2 <- rbinom(n,1,0.5)

 u <- runif(n)
 y <- 1*qexp(u) + (2 + 3*u)*x1 + 5*x2

 # This model is likely to suffer from quantile crossing
 m <- iqr(y ~ x1 + x2, formula.p = ~ slp(p,7))
 diagnose.qc(m)

 # Repeat estimation with remove.qc = TRUE
 m2 <- iqr(y ~ x1 + x2, formula.p = ~ slp(p,7), remove.qc = TRUE)
 diagnose.qc(m2)

 # Use remove.qc = qc.control(trace = TRUE) to see what is going on!
 # You can set a larger 'maxTry', if the algorithm failed to remove 
 # quantile crossing entirely, or a smaller one, if you want to stop
 # the procedure before it becomes 'too expensive' in terms of loss.
 
}
}
\keyword{methods}

