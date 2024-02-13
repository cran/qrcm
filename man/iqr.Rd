\name{iqr}
\alias{iqr}
\title{
Quantile Regression Coefficients Modeling
}
\description{
This function implements Frumento and Bottai's (2016, 2017) and Hsu, Wen, and Chen's (2021) methods for quantile regression coefficients modeling (qrcm). Quantile regression coefficients are described by (flexible) parametric functions of the order of the quantile. Quantile crossing can be eliminated using the method described in Sottile and Frumento (2023).
}
\usage{
iqr(formula, formula.p = ~ slp(p,3), weights, data, s, 
    tol = 1e-6, maxit, remove.qc = FALSE)
}
\arguments{
  \item{formula}{
    a two-sided formula of the form \code{y ~ x1 + x2 + \ldots}: 
    a symbolic description of the quantile regression model. The left side of the formula
    is \code{Surv(time,event)} if the data are right-censored; \code{Surv(time,time2,event)}
    if the data are right-censored and left-truncated (\code{time < time2}, \code{time} can be \kbd{-Inf});
    and \code{Surv(time1, time2, type = "interval2")} for interval-censored data (use \code{time1 = time2} 
    for exact observations, \code{time1 = -Inf} or \code{NA} for left-censored, 
    and \code{time2 = Inf} or \code{NA} for right-censored).
  }
  \item{formula.p}{
    a one-sided formula of the form \code{~ b1(p, \ldots) + b2(p, \ldots) + \ldots}, describing how 
    quantile regression coefficients depend on \kbd{p}, the order of the quantile.
  }
  \item{weights}{
    an optional vector of weights to be used in the fitting process. The weights will always be normalized
    to sum to the sample size. This implies that, for example, using double weights will \emph{not} halve the standard errors.}
  \item{data}{
    an optional data frame, list or environment containing the variables in \code{formula}.
  }
  \item{s}{an optional 0/1 matrix that permits excluding some model coefficients
  (see \sQuote{Examples}).}
  \item{tol}{
    convergence criterion for numerical optimization.
  }
  \item{maxit}{
    maximum number of iterations.
  }
  \item{remove.qc}{
    either a logical value, or a list created with \code{\link{qc.control}}. See \sQuote{Details}.
  }
}
\details{

  Quantile regression permits modeling conditional quantiles of a response variabile,
  given a set of covariates. A linear model is used to describe the conditional
  quantile function:
  \deqn{Q(p | x) = \beta_0(p) + \beta_1(p)x_1 + \beta_2(p)x_2 + \ldots.}{%
  Q(p | x) = \beta0(p) + \beta1(p)*x1 + \beta2(p)*x2 + \ldots.}
  The model coefficients \eqn{\beta(p)} describe the effect of covariates on the \eqn{p}-th
  quantile of the response variable. Usually, one or more
  quantiles  are estimated, corresponding to different values of \eqn{p}.
  
  Assume that each coefficient can be expressed as a parametric function of \eqn{p} of the form:
  \deqn{\beta(p | \theta) = \theta_{0} + \theta_1 b_1(p) + \theta_2 b_2(p) + \ldots}{%
  \beta(p | \theta) = \theta0 + \theta1*b1(p) + \theta2*b2(p) + \ldots}
  where \eqn{b_1(p), b_2(p, \ldots)}{b1(p), b2(p), \ldots} are known functions of \eqn{p}. 
  If \eqn{q} is the dimension of 
  \eqn{x = (1, x_1, x_2, \ldots)}{x = (1, x1, x2, \ldots)} 
  and \eqn{k} is that of  
  \eqn{b(p) = (1, b_1(p), b_2(p), \ldots)}{b(p) = (1, b1(p), b2(p), \ldots)}, 
  the entire conditional quantile function is described by a 
  \eqn{q \times k}{q*k} matrix \eqn{\theta} of model parameters. 
  
  
  Users are required to specify two formulas: \code{formula} describes the regression model,
  while \code{formula.p} identifies the 'basis' \eqn{b(p)}. 
  By default, \code{formula.p = ~ slp(p, k = 3)}, a 3rd-degree shifted 
  Legendre polynomial (see \code{\link{slp}}). Any user-defined function \eqn{b(p, \ldots)}
  can be used, see \sQuote{Examples}.
  
  If no censoring and truncation are present, estimation of \eqn{\theta} is carried out 
  by minimizing an objective function that corresponds 
  to the integral, with respect to \eqn{p}, of the loss function of standard quantile regression. 
  Details are in Frumento and Bottai (2016). If the data are censored or truncated, instead,
  \eqn{\theta} is estimated by solving the estimating equations described in Frumento and Bottai (2017)
  and Hsu, Wen, and Chen (2021).

  The option \code{remove.qc} applies the method described by Sottile and Frumento (2023) 
  to remove quantile crossing. You can either choose \code{remove.qc = TRUE}, or use
  \code{remove.qc = qc.control(...)}, which allows to specify the operational parameters
  of the algorithm. Please read \code{\link{qc.control}} for more details on the method,
  and use \code{\link{diagnose.qc}} to diagnose quantile crossing.
}
\value{
An object of class \dQuote{\code{iqr}}, a list containing the following items:
\item{coefficients}{a matrix of estimated model parameters describing the fitted quantile function.}
\item{converged}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{call}{the matched call.}
\item{obj.function}{if the data are neither censored nor truncated, the value of the minimized loss function; otherwise, a meaningful loss function which, however, is not the objective function of the model (see note 3). The number of model parameter is returned as an attribute.}
\item{mf}{the model frame used.}
\item{PDF, CDF}{the fitted values of the conditional probability density function (\acronym{PDF})
and cumulative distribution function (\acronym{CDF}). See note 1 for details.}
\item{covar}{the estimated covariance matrix.}
\item{s}{the used \sQuote{\kbd{s}} matrix.}
Use \code{\link{summary.iqr}}, \code{\link{plot.iqr}}, and \code{\link{predict.iqr}}
for summary information, plotting, and predictions from the fitted model. The function
\code{\link{test.fit}} can be used for goodness-of-fit assessment.
The generic accessory functions \code{coefficients}, \code{formula}, \code{terms}, \code{model.matrix}, 
\code{vcov} are available to extract information from the fitted model. The special function
\code{\link{diagnose.qc}} can be used to diagnose quantile crossing.
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), 74-84.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics, 73 (4), 1179-1188.

Frumento, P., and Salvati, N. (2021). \emph{Parametric modeling of quantile regression coefficient functions with count data}. Statistical Methods and Applications, 30, 1237-1258.

Hsu, C.Y., Wen, C.C., and Chen, Y.H. (2021). \emph{Quantile function regression analysis for interval censored data, with application to salary survey data}. Japanese Journal of Statistics and Data Science, 4, 999-1018.

Sottile, G., and Frumento, P. (2023). \emph{Parametric estimation of non-crossing quantile functions}. Statistical Modelling, 23 (2), 173-195.

Frumento, P., and Corsini, L. (2024). \emph{Using parametric quantile regression
to investigate determinants of unemployment duration}. Unpublished manuscript.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\note{

  \strong{NOTE 1 (\kbd{PDF}, \kbd{CDF}, quantile crossing, and goodness-of-fit).} 
  By expressing quantile regression coefficients as functions of \eqn{p}, you practically specify 
  a parametric model for the entire conditional distribution. The induced \kbd{CDF} is the value \eqn{p^*}
  such that \eqn{y = Q(p^* | x)}. The corresponding \kbd{PDF} is given by \eqn{1/Q'(p^* | x)}.
  Negative values of \code{PDF} indicate quantile crossing, occurring when the estimated quantile function is not
  monotonically increasing. If negative \code{PDF} values occur for a relatively large proportion of data, 
  the model is probably misspecified or ill-defined.
  If the model is correct, the fitted \code{CDF} should approximately follow a Uniform(0,1) distribution.
  This idea is used to implement a goodness-of-fit test, see \code{\link{test.fit}}.

  \strong{NOTE 2 (model intercept).} 
  The intercept can be excluded from \code{formula}, e.g., 
  \code{iqr(y ~ -1 + x)}. This, however, implies that when \code{x = 0},
  \code{y} is zero at all quantiles. See example 5 in \sQuote{Examples}.
  The intercept can also be removed from \code{formula.p}. 
  This is recommended if the data are bounded. For example, for strictly positive data,
  use \code{iqr(y ~ 1, formula.p = -1 + slp(p,3))} to force the smallest quantile
  to be zero. See example 6 in \sQuote{Examples}.
  
  \strong{NOTE 3 (censoring, truncation, and loss function).} 
  Data are right-censored when, instead of a response variable \eqn{T}, one can only observe
  \eqn{Y = min(T,C)} and \eqn{d = I(T \le C)}. Here, \eqn{C} is a censoring variable
  that is assumed to be conditionally independent of \eqn{T}. Additionally, left truncation occurs if
  \eqn{Y} can only be observed when it exceeds another random variable \eqn{Z}. For example,
  in the prevalent sampling design, subjects with a disease are enrolled; those who died
  before enrollment are not observed. 
  
  Ordinary quantile regression minimizes \eqn{L(\beta(p)) = \sum (p - \omega)(t - x'\beta(p))}
  where \eqn{\omega = I(t \le x'\beta(p))}. Equivalently, it solves its first derivative,
  \eqn{S(\beta(p)) = \sum x(\omega - p)}. The objective function of \command{iqr}
  is simply the integral of \eqn{L(\beta(p | \theta))} with respect to \eqn{p}.
  
  If the data are censored and truncated, \eqn{\omega} is replaced by
  \deqn{\omega^* = \omega.y + (1 - d)\omega.y(p - 1)/S.y - \omega.z - \omega.z(p - 1)/S.z + p}
  where \eqn{\omega.y = I(y \le x'\beta(p))}, \eqn{\omega.z = I(z \le x'\beta(p))}, \eqn{S.y = P(T > y)},
  and \eqn{S.z = P(T > z)}. 
  The above formula can be obtained from equation (7) of Frumento and Bottai, 2017. 
  Replacing \eqn{\omega} with \eqn{\omega^*} in \eqn{L(\beta(p))} is \strong{NOT} equivalent 
  to replacing \eqn{\omega} with \eqn{\omega^*} in \eqn{S(\beta(p))}.
  % Note to myself: I am totally sure of this! The two functions L(omega = omega*) and S(omega = omega*) 
  % are NOT the derivative of each other, and do NOT generate the same solution, not even asymptotically.
  The latter option leads to a much simpler computation, and generates the estimating 
  equation used by \command{iqr}. This means that, if the data are censored or truncated, 
  the \kbd{obj.function} returned by \command{iqr} is \strong{NOT} the objective function being
  minimized, and should not be used to compare models. However, if one of two models has a much larger
  value of the \kbd{obj.function}, this may be a sign of severe misspecification or poor convergence.
  
  If the data are interval-censored, the loss function is obtained as the average between the loss calculated on the lower end of the interval, and that calculated on the upper end. The presence of right- or left-censored observations is handled as described above.
}
\seealso{
\code{\link{summary.iqr}}, \code{\link{plot.iqr}}, \code{\link{predict.iqr}},
for summary, plotting, and prediction, and \code{\link{test.fit.iqr}} for goodness-of-fit assessment;
\code{\link{plf}} and \code{\link{slp}} to define \eqn{b(p)}
to be a piecewise linear function and a shifted Legendre polynomial basis, respectively;
\code{\link{diagnose.qc}} to diagnose quantile crossing.
}
\examples{

  ##### Using simulated data in all examples


  ##### Example 1
  
  n <- 1000
  x <- runif(n)
  y <- rnorm(n, 1 + x, 1 + x)
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with 
    # beta0(p) = beta1(p) = 1 + qnorm(p)
                              
  # fit the true model: b(p) = (1 , qnorm(p))
  m1 <- iqr(y ~ x, formula.p = ~ I(qnorm(p)))
  # the fitted quantile regression coefficient functions are
    # beta0(p) = m1$coef[1,1] + m1$coef[1,2]*qnorm(p)
    # beta1(p) = m1$coef[2,1] + m1$coef[2,2]*qnorm(p)
  \donttest{
  # a basis b(p) = (1, p), i.e., beta(p) is assumed to be a linear function of p
  m2 <- iqr(y ~ x, formula.p = ~ p)

  # a 'rich' basis b(p) = (1, p, p^2, log(p), log(1 - p))
  m3 <- iqr(y ~ x, formula.p = ~ p + I(p^2) + I(log(p)) + I(log(1 - p)))

  # 'slp' creates an orthogonal spline basis using shifted Legendre polynomials
  m4 <- iqr(y ~ x, formula.p = ~ slp(p, k = 3)) # note that this is the default
  
  # 'plf' creates the basis of a piecewise linear function
  m5 <- iqr(y ~ x, formula.p = ~ plf(p, knots = c(0.1,0.9)))

  
  summary(m1)
  summary(m1, p = c(0.25,0.5,0.75))
  test.fit(m1)
  par(mfrow = c(1,2)); plot(m1, ask = FALSE)
  # see the documentation for 'summary.iqr', 'test.fit.iqr', and 'plot.iqr'
  }



  ##### Example 2 ### excluding coefficients
  
  n <- 1000
  x <- runif(n)
  qy <- function(p,x){(1 + qnorm(p)) + (1 + log(p))*x}
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
    # beta0(p) = 1 + qnorm(p) 
    # beta1(p) = 1 + log(p)
  
  y <- qy(runif(n), x) # to generate y, plug uniform p in qy(p,x) 
  iqr(y ~ x, formula.p = ~ I(qnorm(p)) + I(log(p)))

  # I would like to exclude log(p) from beta0(p), and qnorm(p) from beta1(p)
  # I set to 0 the corresponding entries of 's'

  s <- matrix(1,2,3); s[1,3] <- s[2,2] <- 0
  iqr(y ~ x, formula.p = ~ I(qnorm(p)) + I(log(p)), s = s)


  \donttest{

  ##### Example 3 ### excluding coefficients when b(p) is singular
  
  n <- 1000
  x <- runif(n)
  qy <- function(p,x){(1 + log(p) - 2*log(1 - p)) + (1 + log(p/(1 - p)))*x} 
  # true quantile function: Q(p | x) = beta0(p) + beta1(p)*x, with
    # beta0(p) = 1 + log(p) - 2*log(1 - p)
    # beta1(p) = 1 + log(p/(1 - p))

  y <- qy(runif(n), x) # to generate y, plug uniform p in qy(p,x)

  iqr(y ~ x, formula.p = ~ I(log(p)) + I(log(1 - p)) + I(log(p/(1 - p))))
  # log(p/(1 - p)) is dropped due to singularity
  
  # I want beta0(p) to be a function of log(p) and log(1 - p),
  # and beta1(p) to depend on log(p/(1 - p)) alone

  s <- matrix(1,2,4); s[2,2:3] <- 0
  iqr(y ~ x, formula.p = ~ I(log(p)) + I(log(1 - p)) + I(log(p/(1 - p))), s = s)
  # log(p/(1 - p)) is not dropped




  ##### Example 4 ### using slp to test deviations from normality
  
  n <- 1000
  x <- runif(n)
  y <- rnorm(n, 2 + x) 
  # the true model is normal, i.e., b(p) = (1, qnorm(p))
  
  summary(iqr(y ~ x, formula.p = ~ I(qnorm(p)) + slp(p,3))) 
  # if slp(p,3) is not significant, no deviation from normality




  ##### Example 5 ### formula without intercept
  
  n <- 1000
  x <- runif(n)
  y <- runif(n, 0,x) 

  # True quantile function: Q(p | x) = p*x, i.e., beta0(p) = 0, beta1(p) = p
  # When x = 0, all quantiles of y are 0, i.e., the distribution is degenerated
  # To explicitly model this, remove the intercept from 'formula'
  
  iqr(y ~ -1 + x, formula.p = ~ p)
  
  # the true model does not have intercept in b(p) either:

  iqr(y ~ -1 + x, formula.p = ~ -1 + p)




  ##### Example 6 ### no covariates, strictly positive outcome
  
  n <- 1000
  y <- rgamma(n, 3,1) 

  # you know that Q(0) = 0
  # remove intercept from 'formula.p', and use b(p) such that b(0) = 0
  
  summary(iqr(y ~ 1, formula.p = ~ -1 + slp(p,5))) # shifted Legendre polynomials
  summary(iqr(y ~ 1, formula.p = ~ -1 + sin(p*pi/2) + I(qbeta(p,2,4)))) # unusual basis
  summary(iqr(y ~ 1, formula.p = ~ -1 + I(sqrt(p))*I(log(1 - p)))) # you can include interactions




  ##### Example 7 ### revisiting the classical linear model
  
  n <- 1000
  x <- runif(n)
  y <- 2 + 3*x + rnorm(n,0,1) # beta0 = 2, beta1 = 3
  
  iqr(y ~ x, formula.p = ~ I(qnorm(p)), s = matrix(c(1,1,1,0),2))
  # first column of coefficients: (beta0, beta1)
  # top-right coefficient: residual standard deviation
  
  }




  ##### Example 8 ### censored data
  
  n <- 1000
  x <- runif(n,0,5)
	
  u <- runif(n)
  beta0 <- -log(1 - u)
  beta1 <- 0.2*log(1 - u)
  t <- beta0 + beta1*x  # time variable
  c <- rexp(n,2)        # censoring variable
  y <- pmin(t,c)        # observed events
  d <- (t <= c)         # 1 = event, 0 = censored
  
  iqr(Surv(y,d) ~ x, formula.p = ~ I(log(1 - p)))
  
  ##### Example 8 (cont.) ### censored and truncated data

  z <- rexp(n,10)   # truncation variable
  w <- which(y > z) # only observe z,y,d,x when y > z
  z <- z[w]; y <- y[w]; d <- d[w]; x <- x[w]

  iqr(Surv(z,y,d) ~ x, formula.p = ~ I(log(1 - p)))
  
  
  
  
  ##### Example 9 ### interval-censored data
  # (with a very naif data-generating process)
  
  n <- 1000
  x <- runif(n,0,5)
	
  u <- runif(n)
  beta0 <- 10*u + 20*u^2
  beta1 <- 10*u
  t <- beta0 + beta1*x  # time variable
  time1 <- floor(t)     # lower bound
  time2 <- ceiling(t)    # upper bound
  iqr(Surv(time1, time2, type = "interval2") ~ x, formula.p = ~ -1 + p + I(p^2))
  
}
\keyword{models}
\keyword{regression}
