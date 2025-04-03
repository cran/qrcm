\name{qrcm-package}
\alias{qrcm-package}
\docType{package}
\title{
Quantile Regression Coefficients Modeling
}
\description{
This package implements quantile regression coefficient modeling (qrcm), in which the coefficients of a quantile regression model are described by (flexible) parametric functions. The method is described in Frumento and Bottai (2016, 2017); Frumento and Salvati (2021); Frumento, Bottai, and Fernandez-Val (2021); and Hsu, Wen, and Chen (2021). Special functions can be used to diagnose and eliminate quantile crossing (Sottile and Frumento, 2023).
}
\details{
\tabular{ll}{
Package: \tab qrcm\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2025-04-03\cr
License: \tab GPL-2\cr
}
The function \code{\link{iqr}} permits specifying regression models for cross-sectional data, allowing for censored and truncated outcomes. The function \code{\link{iqrL}} can be used to analyze longitudinal data in which the same individuals are observed repeatedly.

Two special functions, \code{\link{slp}} and \code{\link{plf}}, can be used for model building. Auxiliary functions for model summary, prediction, and plotting are provided.
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
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), 74-84.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics, 73 (4), 1179-1188.

Frumento, P., and Salvati, N. (2021). \emph{Parametric modeling of quantile regression coefficient functions with count data}. Statistical Methods and Applications, 30, 1237-1258.

Frumento, P., Bottai, M., and Fernandez-Val, I. (2021). \emph{Parametric modeling of quantile regression coefficient functions with longitudinal data}. Journal of the American Statistical Association, 116 (534), 783-797.

Hsu, C.Y., Wen, C.C., and Chen, Y.H. (2021). \emph{Quantile function regression analysis for interval censored data, with application to salary survey data}. Japanese Journal of Statistics and Data Science, 4, 999-1018.

Sottile, G., and Frumento, P. (2023). \emph{Parametric estimation of non-crossing quantile functions}. Statistical Modelling, 23 (2), 173-195.

Corsini, L., and Frumento, P. (2025). \emph{Using parametric quantile regression to investigate determinants of unemployment duration}. Unpublished manuscript.

}

\keyword{ package }
\examples{
\donttest{
 # iqr(y ~ x) # cross-sectional observations
 # iqr(Surv(time, event) ~ x) # right-censored data
 # iqr(Surv(start, stop, event) ~ x) # right-censored and left-truncated data
 # iqr(Surv(time1, time2, type = "interval") ~ x) # interval-censored data
 # iqrL(y ~ x, id = id) # repeated measures
 
 # diagnose.qc(model) # diagnose quantile crossing
 # Use iqr(..., remove.qc = TRUE) to remove crossing
}
}



