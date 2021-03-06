\name{predict.iqr}
\alias{predict.iqr}
\title{
Prediction After Quantile Regression Coefficients Modeling
}
\description{
Predictions from an object of class \dQuote{\code{iqr}}.
}
\usage{
\method{predict}{iqr}(object, type = c("beta", "CDF", "QF", "sim"), newdata, p, se = TRUE, \ldots)
}
\arguments{
  \item{object}{
 an object of class \dQuote{\code{iqr}}, the result of a call to \code{\link{iqr}}.
}
  \item{type}{
  a character string specifying the type of prediction. See \sQuote{Details}.
}
  \item{newdata}{
  an optional data frame in which to look for variables with which to predict. 
  If omitted, the data are used. For \kbd{type = "CDF"}, it must include the response variable.
 Ignored if \kbd{type = "beta"}.
}
  \item{p}{
  a numeric vector indicating the order(s) of the quantile to predict. Only used if 
	\kbd{type = "beta"} or \kbd{type = "QF"}.
}
  \item{se}{
  logical. If \kbd{TRUE} (the default), standard errors of the prediction will be computed. Only used if \kbd{type = "beta"} or \kbd{type = "QF"}.
}
\item{\ldots}{for future methods.}
}
\details{
Using \code{\link{iqr}}, quantile regression coefficients 
\eqn{\beta(p)} are modeled as parametric functions of \eqn{p}, the order of the quantile.
This implies that the model parameter is \emph{not} \eqn{\beta(p)} itself.
The function \command{predict.iqr} permits computing \eqn{\beta(p)} and other
quantities of interest, as detailed below.
\itemize{
\item if \kbd{type = "beta"} (the default), \eqn{\beta(p)} is returned at 
the supplied value(s) of \kbd{p}. If \kbd{p} is missing, a default \kbd{p = (0.01, ..., 0.99)} is used. 
\item if \kbd{type = "CDF"}, the value of the fitted \acronym{CDF} (cumulative distribution function)
and \acronym{PDF} (probability density function) are computed.
\item if \kbd{type = "QF"}, the fitted values \eqn{x'\beta(p)}, corresponding to the 
conditional quantile function, are computed at the supplied values of \kbd{p}. 
\item if \kbd{type = "sim"}, data are simulated from the fitted model. 
To simulate the data, the fitted conditional quantile function is computed
at randomly generated \kbd{p} following a Uniform(0,1) distribution.
}
}
\value{
\itemize{
\item if \kbd{type = "beta"} a list with one item for each covariate in the model. 
Each element of the list is a data frame with columns (\kbd{p, beta, se, low, up}) reporting \eqn{\beta(p)}, its estimated standard error, and the corresponding 95\% confidence interval. If \kbd{se = FALSE}, the last three columns are not computed.
\item if \kbd{type = "CDF"}, a two-columns data frame \kbd{(CDF,PDF)}.
\item if \kbd{type = "QF"} and \kbd{se = FALSE}, a data frame with one row
for each observation, and one column for each value of \kbd{p}. If \kbd{se = TRUE},
a list of two data frames, \kbd{fit} (predictions) and \kbd{se.fit} (standard errors).
\item if \kbd{type = "sim"}, a vector of simulated data.
}}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\note{
Prediction may generate quantile crossing
if the support of the new covariates values supplied in \code{newdata} 
is different from that of the observed data.
}
\seealso{
\code{\link{iqr}}, for model fitting; \code{\link{summary.iqr}} and \code{\link{plot.iqr}},
for summarizing and plotting \code{iqr} objects.
}
\examples{

  # using simulated data
  
  n <- 1000
  x <- runif(n)
  y <- rlogis(n, 1 + x, 1 + x)
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
    # beta0(p) = beta1(p) = 1 + log(p/(1 - p))
  
  model <- iqr(y ~ x, formula.p = ~ I(log(p)) + I(log(1 - p))) 
  # (fit asymmetric logistic distribution)
  
  
  # predict beta(0.25), beta(0.5), beta(0.75)
  predict(model, type = "beta", p = c(0.25,0.5, 0.75))
  
  # predict the CDF and the PDF at new values of x and y
  predict(model, type = "CDF", newdata = data.frame(x = c(.1,.2,.3), y = c(1,2,3)))
  
  # computes the quantile function at new x, for p = (0.25,0.5,0.75)
  predict(model, type = "QF", p = c(0.25,0.5,0.75), newdata = data.frame(x = c(.1,.2,.3)))

  # simulate data from the fitted model
  ysim <- predict(model, type = "sim") # 'newdata' can be supplied

  # if the model is correct, the distribution of y and that of ysim should be similar
  qy <- quantile(y, prob = seq(.1,.9,.1))
  qsim <- quantile(ysim, prob = seq(.1,.9,.1))
  plot(qy, qsim); abline(0,1)
}
\keyword{methods}
