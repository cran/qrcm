\name{qrcm-package}
\alias{qrcm-package}
\docType{package}
\title{
Quantile Regression Coefficients Modeling
}
\description{
This package implements Frumento and Bottai's (2016, 2017) method for quantile regression
coefficient modeling (qrcm), in which quantile regression coefficients are described
by (flexible) parametric functions of the order of the quantile. This permits modeling
the entire conditional quantile function of a response variable. 
}
\details{
\tabular{ll}{
Package: \tab qrcm\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2017-03-26\cr
License: \tab GPL-2\cr
}
The function \code{\link{iqr}} permits specifying the regression model, allowing for censored and truncated outcomes.
Two special functions, \code{\link{slp}} and \code{\link{plf}}, are provided to facilitate model building.
The auxiliary functions \code{\link{summary.iqr}}, \code{\link{predict.iqr}}, and \code{\link{plot.iqr}}
can be used to extract information from the fitted model, while \code{\link{test.fit}}
is used to assess the model fit.
}
\author{
Paolo Frumento

Maintainer: Paolo Frumento <paolo.frumento@ki.se>
}
\references{
Frumento, P., and Bottai, M. (2016). \emph{Parametric modeling of quantile regression coefficient functions}. Biometrics, 72 (1), pp 74-84, doi: 10.1111/biom.12410.

Frumento, P., and Bottai, M. (2017). \emph{Parametric modeling of quantile regression coefficient functions with censored and truncated data}. Biometrics, doi: 10.1111/biom.12675.

}

\keyword{ package }
\examples{

# use simulated data

n <- 1000
x <- rexp(n)
y <- runif(n, 0, 1 + x)
model <- iqr(y ~ x, formula.p = ~ p + I(p^2))
summary(model)
summary(model, p = c(0.1,0.2,0.3))
predict(model, type = "beta", p = c(0.1,0.2,0.3))
predict(model, type = "CDF", newdata = data.frame(x = c(1,2,3), y = c(0.5,1,2)))
predict(model, type = "QF", p = c(0.1,0.2,0.3), newdata = data.frame(x = c(1,2,3)))
predict(model, type = "sim", newdata = data.frame(x = c(1,2,3)))
par(mfrow = c(1,2)); plot(model, ask = FALSE)
test.fit(model, R = 30)
}



