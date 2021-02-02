\name{plot.iqrL}
\alias{plot.iqrL}
\title{
Plot Quantile Regression Coefficients with Longitudinal Data
}
\description{
Plots quantile regression coefficients 
\eqn{\beta(u)} and \eqn{\gamma(v)},
based on a fitted model of class \dQuote{\code{iqrL}}.
}
\usage{
\method{plot}{iqrL}(x, conf.int = TRUE, polygon = TRUE, which = NULL, ask = TRUE, \ldots)
}
\arguments{
  \item{x}{
an object of class \dQuote{\code{iqrL}}, the result of a call to \code{\link{iqrL}}.
}
  \item{conf.int}{
  logical. If \kbd{TRUE}, asymptotic 95\% confidence intervals are added to the plot.
}
  \item{polygon}{
  logical. If \kbd{TRUE}, confidence intervals are represented by shaded areas via \code{polygon}.
  Otherwise, dashed lines are used.
}
  \item{which}{
  an optional numerical vector indicating which coefficient(s) to plot. If \kbd{which = NULL}, 
  all coefficients are plotted.
}
\item{ask}{logical. If \kbd{which = NULL} and \kbd{ask = TRUE} (the default), you will be asked
interactively which coefficients to plot. Additional options will permit
creating Q-Q plots of \code{u} or \code{v}, which should be independently distributed according to a Uniform(0,1) 
distribution. The option \code{ppplot(u,v)} will generate a P-P plot that compares
the empirical distribution of \code{(u,v)} with its theoretical
value, \code{F(u,v) = uv}, at a discrete grid of points.}
\item{\ldots}{additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd, cex.lab, cex.axis, axes, frame.plot}. 
See \code{\link{par}}.}
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\seealso{
\code{\link{iqrL}} for model fitting; \code{\link{summary.iqrL}} and \code{\link{predict.iqrL}}
for model summary and prediction.
}
\examples{
  
  # using simulated data
  
  n <- 1000 # n. of observations
  n.id <- 100 # n. of clusters
  id <- rep(1:n.id, each = n/n.id) # cluster id

  x1 <- runif(n) # a level-1 covariate
  z1 <- rnorm(n.id) # a level-2 covariate

  V <- runif(n.id) # V_i
  U <- runif(n) # U_it

  alpha <- 2*(V - 1) + z1 # alpha
  y_alpha <- 1 + 2*qnorm(U) + 3*U*x1 # y - alpha
  y <- y_alpha + alpha[id] # observed outcome
  mydata <- data.frame(id = id, y = y, x1 = x1, z1 = z1[id])

  model <- iqrL(fx = y ~ x1, fu = ~ I(qnorm(u)) + u, 
    fz = ~ z1, fv = ~ -1 + I(qnorm(v)), id = id, data = mydata)
  par(mfrow = c(2,2))
  plot(model, ask = FALSE) 
}
\keyword{methods}

