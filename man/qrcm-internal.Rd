\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{
ctiqr(formula, formula.p = ~ slp(p,3), weights, data, s, tol = 1e-5, maxit) % to be added later

check.in(mf, formula.p, type, s)
check.out(theta, S, covar)
start.theta(y, z, d, x, weights, bfun, df, xx, yy, zz, s)
ctiqr.internal(mf, cl, formula.p, tol = 1e-06, maxit, type, s)
iqr.ee(theta, y, z, d, X, weights, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
cov.theta(theta, y, z, d, X, weights, bfun, p.star.y, p.star.z, type, s)
ciqr.ee(theta, y, z, d, X, weights, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
ctiqr.ee(theta, y, z, d, X, weights, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
iobjfun(theta, y, X, weights, bfun, p.star)
iqr.newton(theta, y, z, d, X, weights, bfun, type, tol = 1e-05, maxit = 200, s)

numint(f)
apply_bfun(bfun, p, fun, convert = FALSE)
make.bfun(p, x, type)
p.bisec(theta, y, X, bfun, n.it = 17, convert = FALSE)
slp.basis(k, intercept)
is.slp(f)

chitest(y, w, k = 10, type, Fy)
iqr.waldtest(obj)
km(z,y,d,w, type = c("ciqr", "ctiqr"))
extract.p(model, p, cov = FALSE)
pred.beta(model, p, se = FALSE)
\method{print}{iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{terms}{iqr}(x, \ldots)
\method{model.matrix}{iqr}(object, \ldots)
\method{vcov}{iqr}(object, \ldots)
\method{nobs}{iqr}(object, \ldots)
}
\keyword{internal}

