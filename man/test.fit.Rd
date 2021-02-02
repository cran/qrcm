\name{test.fit}
\alias{test.fit}
\title{
Goodness-of-Fit Test
}
\description{
Generic method for goodness-of-fit test. 
}
\usage{
test.fit(object, \ldots)
}
\arguments{
  \item{object}{
    an object of class \dQuote{\code{iqr}} or \dQuote{\code{iqrL}}.
  }
  \item{\ldots}{
   additional arguments to be supplied to \code{\link{test.fit.iqr}} or \code{\link{test.fit.iqrL}}.
  }
}
\details{
  This function will simply call \code{test.fit.iqr} or \code{test.fit.iqrL} depending on \code{class(object)}.
}
\value{
The test statistic(s) and the associated p-values evaluated with Monte Carlo.
}
\seealso{
\code{\link{test.fit.iqr}}, \code{\link{test.fit.iqrL}}
}
\keyword{methods}
