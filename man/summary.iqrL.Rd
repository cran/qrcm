\name{summary.iqrL}
\alias{summary.iqrL}
\title{
Summary After Quantile Regression Coefficients Modeling With Longitudinal Data
}
\description{
Summary of an object of class \dQuote{\code{iqrL}}.
}
\usage{
\method{summary}{iqrL}(object, p, level, cov = FALSE, \ldots)
}
\arguments{
  \item{object}{
   an object of class \dQuote{\code{iqrL}}, the result of a call to \code{\link{iqrL}}.
  }
  \item{p}{
   an optional vector of quantiles.
  }
  \item{level}{
   a numeric scalar. Use \code{level = 1} to summarize \eqn{\beta(u)}, and 
   \code{level = 2} to summarize \eqn{\gamma(v)}. Ignored if \kbd{p} is missing.
  }
  \item{cov}{
  logical. If \kbd{TRUE}, the covariance matrix of the coefficients or is reported.
  Ignored if \kbd{p} is missing.
}
\item{\ldots}{for future methods.}
}

\value{

If \kbd{p} is supplied, 
a standard summary of the estimated quantile regression coefficients 
is returned for each value of \kbd{p}: if \kbd{level = 1},
a summary of \code{beta(p)}, and if \kbd{level = 2}, a summary of \code{gamma(p)}.
If \kbd{cov = TRUE}, the covariance matrix is also reported.

If \kbd{p} is missing (the default), a list with the following items:
\item{converged}{logical value indicating the convergence status.}
\item{n.it}{the number of iterations.}
\item{n}{the number of observations.}
\item{n.id}{the number of unique \kbd{id}s.}
\item{free.par}{the number of free parameters in the model, excluding fixed effects.}
\item{theta}{the estimate of \eqn{\theta}.}
\item{se.theta}{the estimated standard errors associated with \kbd{theta}.}
\item{phi}{the estimate of \eqn{\phi}.}
\item{se.phi}{the estimated standard errors associated with \kbd{phi}.}
\item{test.row.theta, test.row.phi}{Wald test for the covariates. Each \emph{row} of \code{theta} and \code{phi} is tested for nullity.}
\item{test.col.theta, test.col.phi}{Wald test for the building blocks of the quantile function. Each \emph{column} of \code{theta} and \code{phi} is tested for nullity.}
\item{obj.function}{the minimized loss function.}
\item{call}{the matched call.}
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\seealso{
\code{\link{iqrL}}, for model fitting; \code{\link{predict.iqrL}} and \code{\link{plot.iqrL}},
for predicting and plotting objects of class \dQuote{\code{iqrL}}; \code{\link{test.fit.iqrL}}
for a goodness-of-fit test.
}

\examples{

# using simulated data

n <- 1000 # n. of observations
n.id <- 100 # n. of clusters
id <- rep(1:n.id, each = n/n.id) # cluster id

x <- rexp(n) # a covariate

V <- runif(n.id) # V_i
U <- runif(n) # U_it
y <- 1 + 2*log(U) + 3*x + 0.5*qnorm(V)

# true quantile function: Q(u,v | x) = beta0(u) + beta1(u)*x + gamma0(v), with
   # beta0(u) = 1 + 2*log(u)
   # beta1(u) = 3
   # gamma0(v) = 0.5*qnorm(v)

model <- iqrL(fx = y ~ x, fu = ~ 1 + I(log(u)), fz = ~ 1, fv = ~ -1 + I(qnorm(v)), id = id)
summary(model)
summary(model, level = 1, p = c(0.25, 0.75)) # summary of beta(u) at selected quantiles
summary(model, level = 2, p = c(0.1, 0.9)) # summary of gamma(v) at selected quantiles
}
\keyword{methods}
