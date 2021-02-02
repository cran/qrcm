\name{qrcm-package}
\alias{qrcm-package}
\docType{package}
\title{
Quantile Regression Coefficients Modeling
}
\description{
This package implements Frumento and Bottai's (2016, 2017) method for quantile regression
coefficient modeling (qrcm), in which quantile regression coefficients are described
by (flexible) parametric functions of the order of the quantile. The package includes a generalization to longitudinal data (Frumento et al 2021). Special functions can be used to eliminate quantile crossing (Sottile and Frumento 2021).
}
\details{
\tabular{ll}{
Package: \tab qrcm\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2021-01-29\cr
License: \tab GPL-2\cr
}
The function \code{\link{iqr}} permits specifying regression models for cross-sectional data, 
allowing for censored and truncated outcomes. The function \code{\link{iqrL}} can be used to 
analyze longitudinal data in which the same individuals are observed repeatedly.

Two special functions, \code{\link{slp}} and \code{\link{plf}}, can be used for model building.
Auxiliary functions for model summary, prediction, and plotting are provided.
The generic function \code{\link{test.fit}} is used to assess the model fit.

The function \code{\link{diagnose.qc}} can be applied to
\code{iqr} objects to diagnose quantile crossing, and the option \kbd{remove.qc}
can be used to remove it, using the algorithm described in \code{\link{qc.control}}.
}
\author{
Paolo Frumento

Maintainer: Paolo Frumento <paolo.frumento@unipi.it>
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), pp 74-84, doi: 10.1111/biom.12410.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics, doi: 10.1111/biom.12675.

Frumento, P., Bottai, M., and Fernandez-Val, I. (2021). \emph{Parametric modeling of quantile regression coefficient functions with longitudinal data}. Journal of the American Statistical Association [forthcoming].

Sottile, G., and Frumento, P. (2021). \emph{Parametric estimation of non-crossing quantile functions}. Statistical Modelling [forthcoming].
}

\keyword{ package }
\examples{
\donttest{
 # iqr(y ~ x) # cross-sectional observations
 # iqr(Surv(time, event) ~ x) # censored data
 # iqr(Surv(start, stop, event) ~ x) # censored and truncated data
 # iqrL(y ~ x, id = id) # repeated measures
}
}



