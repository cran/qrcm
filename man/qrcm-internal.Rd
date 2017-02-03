\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{
check.in(mf, formula.p, s)
check.out(theta, S, covar)
start.theta(y,z,d, x, weights, bfun, df, yy, zz, s)
ctiqr.internal(mf, cl, formula.p, tol = 1e-06, maxit, s)
iqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
ciqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
ctiqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, J = TRUE, G, i = FALSE)
cov.theta(theta, y, z, d, X, Xw, weights, bfun, p.star.y, p.star.z, type, s)
iobjfun(theta, y, X, weights, bfun, p.star)
iqr.newton(theta, y,z,d,X,Xw, bfun, s, type, tol, maxit, safeit, eps0)

pmax0(x)
num.fun(dx,fx, op = c("int", "der"))
make.bfun(p, x)
apply_bfun(bfun, p, fun = c("bfun", "b1fun"))
p.bisec(theta, y, X, bfun, n.it = 20)
p.bisec.internal(theta, y,X,bp)
slp.basis(k, intercept)
is.slp(f)


iqr.waldtest(obj)
km(z,y,d,w, type, exclude = NULL)
alpha(obj, mz, mc, k = 98, zcmodel, Tc, Tz)
test.unif.ct(z,y,d,w, type, exclude = 0.05)
findagoodestimator(dat,w)
quickpred(obj, y, type = c("PDF", "SF"))
trans(z,y,d,w,type)

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

