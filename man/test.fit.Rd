\name{test.fit}
\alias{test.fit}
\title{
Goodness-of-Fit Test
}
\description{
Goodness-of-fit test for a model
fitted with \code{\link{iqr}}. The Kolmogorov-Smirnov statistic and the Cramer-Von Mises statistic
are computed. Their distribution under the null hypothesis is estimated
with Monte Carlo (see \sQuote{Details}).
}
\usage{
test.fit(object, R = 100)
}
\arguments{
  \item{object}{
    an object of class \dQuote{\code{iqr}}.
  }
  \item{R}{
   number of Monte Carlo replications.
}
}
\details{
  This function permits assessing goodness of fit by testing the null hypothesis
that the \acronym{CDF} values follow a \eqn{U(0,1)} distribution, indicating that
the model is correctly specified.
Since the \acronym{CDF} values depend on estimated parameters, the distribution of 
the test statistic is not known. To evaluate it, the model is fitted on \kbd{R} simulated datasets 
generated under the null hypothesis.
}
\value{
a matrix with columns \code{statistic} and \code{p.value}, 
reporting the Kolmogorov-Smirnov and Cramer-Von Mises statistic and the associated
p-values evaluated with Monte Carlo.
}
\references{
Frumento, P., and Bottai, M. (2015). \emph{Parametric modeling of quantile regression coefficient functions}.
Biometrics, doi: 10.1111/biom.12410.
}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}
\note{A chi-squared test for the same null hypothesis is performed by
\code{\link{summary.iqr}}. We recommend using \command{test.fit} for a more
accurate inference. In particular, the chi-square test is very unreliable
when the model is bounded or has no intercept. 
}

\examples{

y <- rnorm(1000)
m1 <- iqr(y ~ 1, formula.p = ~ I(qnorm(p))) # correct
m2 <- iqr(y ~ 1, formula.p = ~ p)  # misspecified
test.fit(m1)
test.fit(m2)
}
\keyword{htest}
