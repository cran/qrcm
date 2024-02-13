\name{test.fit.iqr}
\alias{test.fit.iqr}
\title{
Goodness-of-Fit Test
}
\description{
Goodness-of-fit test for a model
fitted with \code{\link{iqr}}. The Kolmogorov-Smirnov statistic and the Cramer-Von Mises statistic
are computed. Their distribution under the null hypothesis is evaluated
with Monte Carlo.
}
\usage{
\method{test.fit}{iqr}(object, R = 100, zcmodel, icmodel, trace = FALSE, \ldots)
}
\arguments{
  \item{object}{
    an object of class \dQuote{\code{iqr}}.
  }
  \item{R}{
   number of Monte Carlo replications. If \kbd{R = 0}, the function only returns the test statistics.
  }
  \item{zcmodel}{a numeric value indicating how to model the joint distribution of censoring 
    (\eqn{C}) and truncation (\eqn{Z}). See \sQuote{Details}.
  }
  \item{icmodel}{a list of operational parameters to simulate interval-censored data. See \sQuote{Details}.
  }
  \item{trace}{logical. If \kbd{TRUE}, the progress will be printed.}
  \item{\ldots}{for future arguments.}
}
\details{
  This function permits assessing goodness of fit by testing the null hypothesis
that the \acronym{CDF} values follow a \eqn{U(0,1)} distribution, indicating that
the model is correctly specified.
Since the fitted \acronym{CDF} values depend on estimated parameters, the distribution of 
the test statistic is not known. To evaluate it, the model is fitted on \kbd{R} simulated datasets 
generated under the null hypothesis. 

The testing procedures are described in details by Frumento and Bottai (2016, 2017) and Frumento and Corsini (2024).

\bold{Right-censored and left-truncated data}. If the data are censored and truncated, \code{object$CDF} is as well a censored and truncated outcome, and its quantiles must be computed by using a suitable version of Kaplan-Meier product-limit estimator. The fitted survival curve is then compared with that of a \eqn{U(0,1)} distribution. 

To run Monte Carlo simulations when data are censored or truncated, it is necessary to estimate
the distribution of the censoring and that of the truncation variable. To this goal,
the function \command{pchreg} from the \pkg{pch} package is used, with default settings.

The joint distribution of the censoring variable (\eqn{C}) and the truncation variable (\eqn{Z}) 
can be specified in two ways:
\itemize{
\item If \kbd{zcmodel = 1}, it is assumed that \eqn{C = Z + U},
where \eqn{U} is a positive variable and is independent of \eqn{Z}, given covariates. This is the most
common situation, and is verified when censoring occurs at the end of the follow-up. Under this scenario,
\eqn{C} and \eqn{Z} are correlated with \eqn{P(C > Z) = 1}. 
\item If \kbd{zcmodel = 2}, it is assumed that \eqn{C} and \eqn{Z} are conditionally independent.
This situation is more plausible when all censoring is due to drop-out.
}


\bold{Interval-censored data}. 

If the data are interval-censored, \code{object$CDF} is composed of two columns, \code{left} and \code{right}. A nonparametric estimator is applied to the interval-censored pair \code{(left, right)} using the \pkg{icenReg} R package. The fitted quantiles are then compared with those of a \eqn{U(0,1)} distribution.

To simulate interval-censored data, additional information is required about the censoring mechanism. This testing procedure assumes that interval censoring occurs because each individual is only examined at discrete time points, say \code{t[1], t[2], t[3],}... If this is not the mechanism that generated your data, you should not use this function.

In the ideal situation, one can use \code{t[1], t[2], t[3],}... to estimate the distribution of the time between visits, \code{t[j + 1] - t[j]}. If, however, one only knows \code{time1} and \code{time2}, the two endpoints of the interval, things are more complicated. The empirical distribution of \code{time2 - time1} is NOT a good estimator of the distribution of \code{t[j + 1] - t[j]}, because the events are likely contained in longer intervals, a fact that obviously generates selection bias. There are two common situations: either \code{t[j + 1] - t[j]} is a constant (e.g., one month), or it is random. If \code{t[j + 1] - t[j]} is random and has an Exponential distribution with scale \code{lambda}, then \code{time2 - time1} has a \code{Gamma(shape = 2, scale = lambda)} distribution. This is due to the property of memoryless of the Exponential distribution, and may only be an approximation if there is a floor effect (i.e., if \code{lambda} is larger than the low quantiles of the time-to-event).

The \code{icmodel} argument must be a list with four elements, \code{model}, \code{lambda} (optional), \code{t0}, and \code{logscale}:
\itemize{
\item \code{model}. A character string, either \code{'constant'} or \code{'exponential'}.

\item \code{lambda}. If \code{model = 'constant'}, \code{lambda} will be interpreted as a constant time between visits. If \code{model = 'exponential'}, instead, it will be interpreted as the \emph{mean} (not the rate) of the Exponential distribution that is assumed to describe the time between visits.

If you either \emph{know} \code{lambda}, or you can estimate it by using additional information (e.g., individual data on \emph{all} visit times \code{t[1], t[2], t[3], }...), you can supply a scalar value, that will be used for all individuals, or a vector, allowing \code{lambda} to differ across individuals.

If, instead, \code{lambda} is not supplied or is \code{NULL}, the algorithm proceeds as follows. If \code{model = 'constant'}, the time between visits is assumed to be constant and equal to \code{lambda = mean(time2 - time1)}. If \code{model = 'exponential'}, times between visits are generated from an Exponential distribution in which the mean, \code{lambda}, is allowed to depend on covariates according to a log-linear model, and is estimated by fitting a Gamma model on \code{time2 - time1} as described earlier. 

 

\item \code{t0}. If \code{t0 = 0}, data will be simulated assuming that the first visit occurs at time = 0 (the \dQuote{onset}), i.e., when the individual enters the risk set. This mechanism cannot generate left censoring. If \code{t0 = 1}, instead, the first visit occurs \emph{after} time zero. This mechanism generates left censoring whenever the event occurs before the first visit. Finally, if \code{t0 = -1}, visits start \emph{before} time 0. Under this scenario, it is assumed that not only the time at the event, but also the time at onset is interval-censored. If the event occurs in the interval \code{(time1, time2)}, and the onset is in \code{(t01, t02)}, then the total duration is in the interval \code{(time1 - t02, time2 - t01)}.

\item \code{logscale}. Logical: is the response variable on the log scale? If this is the case, the Monte Carlo procedure will act accordingly. Note that \code{lambda} will always be assumed to describe the time between visits on the \emph{natural} scale.
}

The mechanism described above can automatically account for the presence of left censoring.
In order to simulate right-censored observations (if present in the data), the distribution of the censoring variable is estimated with the function \command{pchreg} from the \pkg{pch} package.


}
\value{
a matrix with columns \code{statistic} and \code{p.value}, 
reporting the Kolmogorov-Smirnov and Cramer-Von Mises statistic and the associated
p-values evaluated with Monte Carlo.
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), pp 74-84, doi: 10.1111/biom.12410.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics, doi: 10.1111/biom.12675.

Frumento, P., and Corsini, L. (2024). \emph{Using parametric quantile regression
to investigate determinants of unemployment duration}. Unpublished manuscript.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}

\examples{
y <- rnorm(1000)
m1 <- iqr(y ~ 1, formula.p = ~ I(qnorm(p))) # correct
m2 <- iqr(y ~ 1, formula.p = ~ p)  # misspecified
\donttest{
test.fit(m1)
test.fit(m2)
}
}
\keyword{htest}
