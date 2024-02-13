\name{iqrL}
\alias{iqrL}
\title{
Quantile Regression Coefficients Modeling with Longitudinal Data
}
\description{
This function implements Frumento et al's (2021) method for quantile regression
coefficients modeling with longitudinal data. 
}
\usage{
iqrL(fx, fu = ~ slp(u,3), fz = ~ 1, fv = ~ -1 + I(qnorm(v)), 
   id, weights, s.theta, s.phi, data, tol = 1e-5, maxit)
}
\arguments{
  \item{fx, fu, fz, fv}{formulas that describe the model (see \sQuote{Details}).}
  \item{id}{a vector of cluster identifiers.}
  \item{weights}{an optional vector of weights to be used in the fitting process.}
  \item{s.theta, s.phi}{optional 0/1 matrices that permit excluding some model coefficients.}
  \item{data}{an optional data frame, list or environment containing the variables in \code{fx} and \code{fz}.}
  \item{tol}{convergence criterion for numerical optimization.}
  \item{maxit}{maximum number of iterations. If missing, a default is computed.}
}
\details{
 
  New users are recommended to read Frumento and Bottai's (2016) paper for details
  on notation and modeling, and to have some familiarity with the 
  \code{\link{iqr}} command, of which \command{iqrL} is a natural expansion.

  The following data-generating process is assumed:
    \deqn{Y_{it} = x_{it}\beta(U_{it}) + z_i\gamma(V_i)}{%
    Y[it] = x[it]*\beta(U[it]) + z[i]*\gamma(V[i])
    }
  where \eqn{x_{it}}{x[it]} are level-1 covariates, \eqn{z_i}{z[i]} are level-2 covariates,
  and \eqn{(U_{it}, V_i)}{(U[it], V[i])} are independent \eqn{U(0,1)} random variables.
  This model implies that \eqn{\alpha_i = z_i\gamma(V_i)}{\alpha[i] = z[i]\gamma(V[i])} are cluster-level
  effects with quantile function \eqn{z_i\gamma(v)}{z[i]\gamma(v)}, while \eqn{x_{it}\beta(u)}{x[it]\beta(u)}
  is the quantile function of \eqn{Y_{it} - \alpha_i}{Y[it] - \alpha[i]}.
  
  Both \eqn{\beta(u)} and \eqn{\gamma(v)}
  are modeled parametrically, using a linear combination of known \dQuote{basis}
  functions \eqn{b(u)} and \eqn{c(v)} such that
    \deqn{\beta(u) = \beta(u | \theta) = \theta b(u),}
    \deqn{\gamma(u) = \gamma(u | \phi) = \phi c(v),}
  where \eqn{\theta} and \eqn{\phi} are matrices of model parameters.

  Model specification is implemented as follows.
  \itemize{
    \item \kbd{fx} is a two-sided formula of the form \kbd{y ~ x}.
    \item \kbd{fu} is a one-sided formula that describes \eqn{b(u)}.
    \item \kbd{fz} is a one-sided formula of the form \kbd{~ z}.
    \item \kbd{fv} is a one-sided formula that describes \eqn{c(v)}.
  }  

  By default, \kbd{fu = ~ slp(u,3)}, a shifted Legendre's polynomial (see \code{\link{slp}}),
  and the distribution of \eqn{\alpha_i}{\alpha[i]} is assumed to be Normal (\kbd{fv = ~ -1 + I(qnorm(v))})
  and to not depend on covariates (\kbd{fz = ~ 1}).

  Restrictions on \eqn{\theta} and \eqn{\phi} are imposed by setting to zero the corresponding elements
  of \kbd{s.theta} and \kbd{s.phi}. 
}
\value{
An object of class \dQuote{\code{iqrL}}, a list containing the following items:
\item{theta, phi}{estimates of \eqn{\theta} and \eqn{\phi}.}
\item{obj.function}{the value of the minimized loss function, and, separately, the level-1 and the level-2 loss. The number of model parameters (excluding the individual effects) is returned as an attribute.}
\item{call}{the matched call.}
\item{converged}{logical. The convergence status.}
\item{n.it}{the number of iterations.}
\item{covar.theta, covar.phi}{the estimated covariance matrices.}
\item{mf.theta, mf.phi}{the model frames used to fit \eqn{\theta} and \eqn{\phi}, respectively. Note that \code{mf.theta} is sorted by increasing \kbd{id} and, within each \kbd{id}, by increasing values of the response variable \kbd{y}, while \code{mf.phi} is sorted by increasing \kbd{id}.}
\item{s.theta, s.phi}{the used \sQuote{\kbd{s.theta}} and \sQuote{\kbd{s.phi}} matrices.}
\item{fit}{
  a data.frame with the following variables:
  \itemize{
    \item \kbd{id} the cluster identifier.
    \item \kbd{y} the response variable.
    \item \kbd{alpha} the estimated individual effects.
    \item \kbd{y_alpha} = \kbd{y - alpha[id]}, the estimated responses purged of the individual effects.
    \item \kbd{v} estimates of \eqn{V_i}{V[i]}.
    \item \kbd{u} estimates of \eqn{U_{it}}{U[it]}.
  }
  Observations are sorted by increasing \kbd{id} and, within each \kbd{id}, by increasing \kbd{y}.
}

Use \code{\link{summary.iqrL}}, \code{\link{plot.iqrL}}, and \code{\link{predict.iqrL}}
for summary information, plotting, and predictions from the fitted model. The function
\code{\link{test.fit.iqrL}} can be used for goodness-of-fit assessment.
The generic accessory functions \code{coefficients}, \code{formula}, \code{terms}, \code{model.matrix}, \code{vcov} are available to extract information from the fitted model.
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), 74-84.

Frumento, P., Bottai, M., and Fernandez-Val, I. (2021). \emph{Parametric modeling of quantile regression coefficient functions with longitudinal data}. Journal of the American Statistical Association, 116 (534), 783-797.
}

\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}

\seealso{
\code{\link{summary.iqrL}}, \code{\link{plot.iqrL}}, \code{\link{predict.iqrL}},
for summary, plotting, and prediction, and \code{\link{test.fit.iqrL}} for goodness-of-fit assessment.
\code{\link{plf}} and \code{\link{slp}} to define \eqn{b(u)} or \eqn{c(v)}
to be piecewise linear functions and shifted Legendre polynomials, respectively.
}
\examples{

  ##### Also see ?iqr for a tutorial on modeling
  ##### Using simulated data in all examples
  

  ##### Example 1

  n <- 1000 # n. of observations
  n.id <- 100 # n. of clusters
  id <- rep(1:n.id, each = n/n.id) # cluster id

  x1 <- runif(n) # a level-1 covariate
  z1 <- rbinom(n,1,0.5)[id] # a level-2 covariate

  V <- runif(n.id) # V_i
  U <- runif(n) # U_it

  alpha <- (0.5 + z1)*qnorm(V) # or alpha = rnorm(n.id, 0, 0.5 + z1)
  y_alpha <- qexp(U) + 3*x1 # or y_alpha = 3*x1 + rexp(n)
  y <- y_alpha + alpha[id] # observed outcome
  mydata <- data.frame(id = id, y = y, x1 = x1, z1 = z1[id])

  # true quantile function: beta0(u) + beta1(u)*x1 + gamma0(v) + gamma1(v)*z1
  # beta0(u) = qexp(u)
  # beta1(u) = 3
  # gamma0(v) = 0.5*qnorm(v)
  # gamma1(v) = qnorm(v)


  ##### Example 1 (cont.) fitting the model

  model1 <- iqrL(fx = y ~ x1, fu = ~ I(qexp(u)), fz = ~ z1, fv = ~ -1 + I(qnorm(v)), 
    id = id, data = mydata)
  summary(model1) # theta, phi
  summary(model1, level = 1, p = c(0.1,0.9)) # beta
  summary(model1, level = 2, p = c(0.1,0.9)) # gamma
  par(mfrow = c(2,2)); plot(model1, ask = FALSE)


  ##### Example 1 (cont.) - excluding coefficients

  s.theta <- rbind(0:1,1:0) # beta0(u) has no intercept, and beta1(u) does not depend on u.
  model2 <- iqrL(fx = y ~ x1, fu = ~ I(qexp(u)), fz = ~ z1, fv = ~ -1 + I(qnorm(v)), 
    id = id, s.theta = s.theta, data = mydata)
  summary(model2)
  \donttest{test.fit(model2)} # testing goodness-of-fit


  ##### Example 1 (cont.) - flexible modeling using slp for lev. 1, asymm. logistic for lev. 2
  \donttest{
  model3 <- iqrL(fx = y ~ x1, fu = ~ slp(u,3), 
    fz = ~ z1, fv = ~ -1 + I(log(2*v)) + I(-log(2*(1 - v))), 
    id = id, data = mydata)
  par(mfrow = c(2,2)); plot(model3, ask = FALSE)


  ##### Example 2 - revisiting the classical linear random-effects model

  n <- 1000 # n. of observations
  n.id <- 100 # n. of clusters
  id <- rep(1:n.id, each = n/n.id) # id

  x1 <- runif(n,0,5)
  E <- rnorm(n) # level-1 error
  W <- rnorm(n.id, 0, 0.5) # level-2 error
  y <- 2 + 3*x1 + E + W[id] # linear random-intercept model

  s.theta <- rbind(1, 1:0)
  linmod <- iqrL(fx = y ~ x1, fu = ~ I(qnorm(u)), id = id, s.theta = s.theta)
  summary(linmod)
  }
}
\keyword{models}
\keyword{regression}
