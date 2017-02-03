\name{test.fit}
\alias{test.fit}
\title{
Goodness-of-Fit Test
}
\description{
Goodness-of-fit test for a model
fitted with \code{\link{iqr}}. The Kolmogorov-Smirnov statistic and the Cramer-Von Mises statistic
are computed. Their distribution under the null hypothesis is estimated
with Monte Carlo.
}
\usage{
test.fit(object, R = 100, zcmodel = 1, trace = FALSE)
}
\arguments{
  \item{object}{
    an object of class \dQuote{\code{iqr}}.
  }
  \item{R}{
   number of Monte Carlo replications.
  }
  \item{zcmodel}{a numeric value indicating how to model the joint distribution of censoring 
    (\eqn{C}) and truncation (\eqn{Z}). Only used when data are censored and truncated. See \sQuote{Details}.
  }
  \item{trace}{logical. If \kbd{TRUE}, the progress be printed.
  }
}
\details{
  This function permits assessing goodness of fit by testing the null hypothesis
that the \acronym{CDF} values follow a \eqn{U(0,1)} distribution, indicating that
the model is correctly specified.
Since the \acronym{CDF} values depend on estimated parameters, the distribution of 
the test statistic is not known. To evaluate it, the model is fitted on \kbd{R} simulated datasets 
generated under the null hypothesis. 

If the data are censored and truncated, \code{object$CDF} is as well a censored and truncated outcome,
and its quantiles are estimated using a null Cox model. The fitted survival curve is then compared with
a \eqn{U(0,1)}. 

To run Monte Carlo simulations when data are censored or truncated,
the distribution of the censoring and that of the truncation variable must be estimated:
the function \command{pchreg} from the \pkg{pch} package is used, with its option \code{splinex = splinex()}.

The joint distribution of the censoring variable (\eqn{C}) and the truncation variable (\eqn{Z}) 
can be specified in two ways:
\itemize{
\item If \kbd{zcmodel = 1} (the default), it is assumed that \eqn{C = Z + U},
where \eqn{U} is a positive variable and is independent of \eqn{Z}, given covariates. This is the most
common situation, and is verified when censoring occurs at the end of the follow-up. Under this scenario,
\eqn{C} and \eqn{Z} are correlated with \eqn{P(C > Z) = 1}. 
\item If \kbd{zcmodel = 2}, it is assumed that \eqn{C} and \eqn{Z} are conditionally independent.
This situation is more plausible when all censoring is due to drop-out.
}
The testing procedure is described in details by Frumento and Bottai (2015, 2016).
}
\value{
a matrix with columns \code{statistic} and \code{p.value}, 
reporting the Kolmogorov-Smirnov and Cramer-Von Mises statistic and the associated
p-values evaluated with Monte Carlo.
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72(1), 74-84, doi: 10.1111/biom.12410.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics (in press), doi:10.1111/biom.12675.

}
\author{
Paolo Frumento \email{paolo.frumento@ki.se}
}

\examples{
y <- rnorm(1000)
m1 <- iqr(y ~ 1, formula.p = ~ I(qnorm(p))) # correct
m2 <- iqr(y ~ 1, formula.p = ~ p)  # misspecified
test.fit(m1)
test.fit(m2)
}
\keyword{htest}
