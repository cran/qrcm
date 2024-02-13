\name{test.fit.iqrL}
\alias{test.fit.iqrL}
\title{
Goodness-of-Fit Test
}
\description{
Goodness-of-fit test for a model
fitted with \code{\link{iqrL}}. The Kolmogorov-Smirnov statistic is computed and its
distribution under the null hypothesis is evaluated with Monte Carlo.
}
\usage{
\method{test.fit}{iqrL}(object, R = 100, trace = FALSE, \ldots)
}
\arguments{
  \item{object}{
    an object of class \dQuote{\code{iqrL}}.
  }
  \item{R}{
   number of Monte Carlo replications. If \kbd{R = 0}, the function only returns the test statistic.
  }
  \item{trace}{logical. If \kbd{TRUE}, the progress will be printed.}
  \item{\ldots}{for future arguments.}
}
\details{
  This function permits assessing goodness of fit by testing the null hypothesis
that the estimated \code{(u,v)} values are independent uniform variables. 
To evaluate the distribution of the test statistic under the true model, a Monte Carlo
method is used (Frumento et al, 2021).
}
\value{
a vector with entries \code{statistic} and \code{p.value}, 
reporting the Kolmogorov-Smirnov statistic (evaluated on a grid)
and the associated p-value.
}
\references{
Frumento, P., Bottai, M., and Fernandez-Val, I. (2021). \emph{Parametric modeling of quantile regression coefficient functions with longitudinal data}. Journal of the American Statistical Association, 116 (534), 783-797.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}

\examples{
id <- rep(1:50, each = 10)
y <- rnorm(500) + rnorm(50)[id]
m1 <- iqrL(fx = y ~ 1, fu = ~ I(qnorm(u)), id = id) # correct
m2 <- iqrL(fx = y ~ 1, fu = ~ u, id = id)  # misspecified
\donttest{
test.fit(m1, R = 20)
test.fit(m2, R = 20)
}
# Warning: this procedure may be time-consuming.
}
\keyword{htest}
