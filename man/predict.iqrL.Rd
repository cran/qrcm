\name{predict.iqrL}
\alias{predict.iqrL}
\title{
Prediction After Quantile Regression Coefficients Modeling with Longitudinal Data
}
\description{
Predictions from an object of class \dQuote{\code{iqrL}}.
}
\usage{
\method{predict}{iqrL}(object, level, type = c("coef", "CDF", "QF", "sim"), newdata, p, se = FALSE, \ldots)
}
\arguments{
  \item{object}{
   an object of class \dQuote{\code{iqrL}}, the result of a call to \code{\link{iqrL}}.
  }
  \item{level}{
   a numeric scalar. Use \code{level = 1} to predict \eqn{y_{it} - \alpha_i}{y[it] - \alpha[i]}, and 
   \code{level = 2} to predict \eqn{\alpha_i}{\alpha[i]} (see \code{\link{iqrL}} for the notation).
  }
  \item{type}{
  a character string specifying the type of prediction. See \sQuote{Details}.
  }
  \item{newdata}{
  an optional data frame in which to look for variables with which to predict (ignored if \kbd{type = "coef"}).
  For \kbd{type = "CDF"}, \code{newdata} must include a response variable named \sQuote{\kbd{y_alpha}}, 
  if \kbd{level = 1}, and \sQuote{\kbd{alpha}} if \kbd{level = 2}. If \code{newdata} is omitted, 
  the observed data will be used, and \kbd{y_alpha} and \kbd{alpha} will be taken from
  \kbd{object$fit}. 
  
  }
  \item{p}{
  a numeric vector indicating the order(s) of the quantile to predict. Only used if 
	\kbd{type = "coef"} or \kbd{type = "QF"}.
  }
  \item{se}{
  logical. If \kbd{TRUE} (the default), standard errors of the prediction will be computed. Only used if \kbd{type = "coef"} or \kbd{type = "QF"}.
  }
  \item{\ldots}{for future methods.}
}
\details{
\itemize{
\item if \kbd{type = "coef"} (the default), quantile regression coefficients are returned:
  if \kbd{level = 1}, \eqn{\beta(p)}; and if \kbd{level = 2}, \eqn{\gamma(p)}. 
  If \kbd{p} is missing, a default \kbd{p = (0.01, ..., 0.99)} is used. 
\item if \kbd{type = "CDF"}, the value of the fitted \acronym{CDF} (cumulative distribution function)
and \acronym{PDF} (probability density function) are computed. If \kbd{level = 1},
these refer to the distribution of \eqn{Y_{it} - \alpha_i = x_{it}\beta(U_{it})}{%
Y[it] - \alpha[i] = x[it]\beta(U[it])},
and the \acronym{CDF} is an estimate of \eqn{U_{it}}{U[it]}. If \kbd{level = 2},
they refer to the distribution of \eqn{\alpha_i = z_i\gamma(V_i)}{\alpha[i] = z[i]\gamma(V[i])},
and the \acronym{CDF} is an estimate of \eqn{V_i}{V[i]}.
\item if \kbd{type = "QF"}, the fitted values \eqn{x\beta(p)} (if \kbd{level = 1}),
or \eqn{z\gamma(p)} (if \kbd{level = 2}). 
\item if \kbd{type = "sim"}, data are simulated from the fitted model. 
If \kbd{level = 1}, simulated values are from the distribution of \eqn{Y_{it} - \alpha_i}{Y[it] - \alpha[i]},
while if \kbd{level = 2}, they are from the distribution of \eqn{\alpha_i}{\alpha[i]}.
}
}
\value{
\itemize{
\item if \kbd{type = "coef"} a list with one item for each covariate. 
Each element of the list is a data frame with columns (\kbd{u, beta, se, low, up}),
if \kbd{level = 1}, and (\kbd{v, gamma, se, low, up}), if \kbd{level = 2}.
If \kbd{se = FALSE}, the last three columns are not computed.
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
If no \kbd{newdata} are supplied, the observed data are used and predictions are ordered as follows:
\itemize{
  \item if \kbd{level = 1}, by increasing \kbd{id} and, within each \kbd{id},
  by increasing values of the response variable \kbd{y}. Rownames will indicate
  the position in the original data frame.
  \item if \kbd{level = 2}, by increasing \kbd{id}.
}
}
\seealso{
\code{\link{iqrL}}, for model fitting; \code{\link{summary.iqrL}} and \code{\link{plot.iqrL}},
for summarizing and plotting \code{iqrL} objects.
}
\examples{

  # using simulated data
  
  n <- 1000 # n. of observations
  n.id <- 100 # n. of clusters
  id <- rep(1:n.id, each = n/n.id) # cluster id

  x1 <- runif(n) # a level-1 covariate
  z1 <- rbinom(n.id,1,0.5) # a level-2 covariate

  V <- runif(n.id) # V_i
  U <- runif(n) # U_it

  alpha <- qlogis(V)*(0.5 + z1) # alpha
  y_alpha <- 1 + 2*qexp(U) + 3*x1 # y - alpha
  y <- y_alpha + alpha[id] # observed outcome
  mydata <- data.frame(id = id, y = y, x1 = x1, z1 = z1[id])

  # true model: Y_it = beta0(U_it) + beta1(U_it)*x1 + gamma0(V_i) + gamma1(V_i)*z1
  # beta0(u) = 1 + 2*pexp(u)
  # beta1(u) = 3
  # gamma0(v) = 0.5*qlogis(v)
  # gamma1(v) = qlogis(V)
  
  model <- iqrL(fx = y ~ x1, fu = ~ I(qexp(u)), fz = ~ z1, fv = ~ -1 + I(qlogis(v)), 
    id = id, data = mydata) 
  
  # predict beta(0.25), beta(0.5), beta(0.75)
  predict(model, level = 1, type = "coef", p = c(0.25,0.5,0.75))
 
  # predict gamma(0.1), gamma(0.9)
  predict(model, level = 2, type = "coef", p = c(0.1,0.9))

  # predict the CDF (u) and the PDF of (y - alpha), at new values of x1
  predict(model, level = 1, type = "CDF", 
    newdata = data.frame(x1 = c(.1,.2,.3), y_alpha = c(1,2,3)))
 
  # predict the CDF (v) and the PDF of alpha, at new values of z1
  predict(model, level = 2, type = "CDF", 
    newdata = data.frame(z1 = c(0,1), alpha = c(-1,1)))
 
  # computes the quantile function of (y - alpha) at new x1, for u = (0.25,0.5,0.75)
  predict(model, level = 1, type = "QF", p = c(0.25,0.5,0.75), 
    newdata = data.frame(x1 = c(.1,.2,.3)))

  # computes the quantile function of alpha at new z1, for v = (0.25,0.5,0.75)
  predict(model, level = 2, type = "QF", p = c(0.25,0.5,0.75), 
    newdata = data.frame(z1 = c(.1,.2,.3)))


  # simulate data from the fitted model
  y_alpha_sim <- predict(model, level = 1, type = "sim")
  alpha_sim <- predict(model, level = 2, type = "sim")
  y_sim = y_alpha_sim + alpha_sim[id]
}
\keyword{methods}
