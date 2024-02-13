\name{diagnose.qc}
\alias{diagnose.qc}
\title{
Diagnose Quantile Crossing
}
\description{
Diagnose quantile crossing in a model estimated with \code{\link{iqr}}.
}
\usage{
diagnose.qc(obj)
}
\arguments{
  \item{obj}{
  an object created with \code{\link{iqr}}.
}
}
\details{
The function determines if quantile crossing occurs in your fitted model, and provides a number of diagnostic tools.

\emph{Local} quantile crossing is defined by \code{obj$PDF < 0}, and is obtained when the quantile function, say \eqn{Q(p | x)}, has negative first derivatives at the values of \eqn{p} that correspond to the observed data.
\emph{Global} quantile crossing occurs when the conditional quantile function has negative first derivatives at \emph{some} values of \eqn{p}. To assess global crossing, a grid of approximately 1000 quantiles is used. Note that local crossing is a special case of global crossing.

The function will assess local and global crossing, and return a summary \kbd{pcross} of the quantiles at which \emph{global} crossing occurs. It is important to understand that crossing at extremely low or high quantiles is very common, but may be considered irrelevant in practice. For example, if \emph{all} observations have crossing quantiles, implying that global crossing is 100\%, but crossing only occurs at quantile above  0.999, the fitted model can be safely used for prediction. Very frequently, crossing occurs at extreme quantiles that do not correspond to any observation in the data. 

This command will also compute a \kbd{crossIndex}, that represents the average length, across observations, of the sub-intervals \eqn{p*} such that \eqn{Q'(p* | x) < 0}. For example, if \eqn{Q'(p | x) < 0} in the interval \eqn{p* = (0.3,0.5)}, the contribution to the \kbd{crossIndex} is \kbd{0.5 - 0.3 = 0.2}. If crossing is detected at a single quantile, the interval is assumed to have length \kbd{1e-6}. In principle, the \kbd{crossIndex} is always between \kbd{0} (no quantile crossing) and \kbd{1} (\emph{all} observations crossing at \emph{all} quantiles, which is clearly impossible). In practice, values of \kbd{crossIndex} greater than 0.05 are relatively rare.

}
\value{
A list with the following items:
\item{qc}{a data frame with two columns \kbd{(qc.local, qc.global)} containing logical indicators of local and global quantile crossing for each observation in the data.}
\item{qc.local, qc.global}{the absolute number of observations for which local/global quantile crossing was detected.}
\item{pcross}{a frequency table of the values of \eqn{p} at which global quantile crossing was detected.}
\item{crossIndex}{the estimated index of crossing described above.}

If no quantile crossing is detected, \code{pcross = NULL}, and \code{crossIndex = 0}. 
}
\references{
Sottile, G., and Frumento, P. (2023). \emph{Parametric estimation of non-crossing quantile functions}. Statistical Modelling, 23(2), 173-195.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\seealso{
\code{\link{iqr}}, \code{\link{qc.control}}.
}
\examples{
 # Using simulated data

 n <- 1000
 x1 <- runif(n,0,3)
 x2 <- rbinom(n,1,0.5)

 u <- runif(n)
 y <- 1*qexp(u) + (2 + 3*u)*x1 + 5*x2
 m <- iqr(y ~ x1 + x2, formula.p = ~ slp(p,7))
 diagnose.qc(m)
}
\keyword{methods}

