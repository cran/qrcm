\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{

check.in.iqr(mf, formula.p, s)
check.in0.iqrL(mf1, mf2)
check.in.iqrL(mf, y,w,formula.p, s)
check.out(theta, S, covar)
start.iqr(y,z,d, x, weights, bfun, df, yy, zz, s)
start.iqrL(y,x,z,id, w1,w2, bfun1,bfun2, s.theta,s.phi, S1,S2)

ctiqr.internal(mf, cl, formula.p, tol = 1e-06, maxit, s, remove.qc)
iqrL.internal(mf1, mf2, cl, fu,fv, s.theta, s.phi, tol = 1e-5, maxit)
iqrL.fit(theta,phi, y,alpha, x,xw,z,zw, id, w1,w2, bfun1,bfun2, s.theta, s.phi, 
	maxit.theta, safeit.theta, maxit.phi, safeit.phi, eps, tol, maxit)

iqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, 
  J = TRUE, G, i = FALSE, lambda = 0)
ciqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, 
  J = TRUE, G, i = FALSE, lambda = 0)
ctiqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, 
  J = TRUE, G, i = FALSE, lambda = 0)
iobjfun(theta, y, X, weights, bfun, p.star)
iobjfun.ct(theta, z,y,d,X,weights, bfun, py, pz, type)
cov.fun.iqr(theta, y, z, d, X, Xw, weights, bfun, p.star.y, p.star.z, type, s)

iqrL.ee(par, x,xw, bfun, p, g = TRUE, H = TRUE, i = FALSE)
cov.fun.iqrL(fit, x,xw,z,zw, id, w1,w2, bfun1,bfun2, s.theta, s.phi)

iqr.newton(theta, y,z,d,X,Xw, bfun, s, type, tol, maxit, safeit, eps0, lambda = 0)
divide.et.impera(fit, V, bfun, s, type, tol, maxit, safeit, eps0, lambda = 0)
iqrL.newton(par, y,x,xw, bfun, s, tol, maxit, safeit, eps)

pmax0(x)
maxind(A)
num.fun(dx,fx, op = c("int", "der"))
make.bfun(p, x)
apply_bfun(bfun, p, fun = c("bfun", "b1fun"))
p.bisec(theta, y, X, bfun, n.it = 20)
p.bisec.internal(theta, y,X,bp)
slp.basis(k, intercept)
is.slp(f)
safesolve(A,B,lambda)

iqr.waldtest(obj)
extract.p(model, p, cov = FALSE)
pred.beta(model, p, se = FALSE)


km(z,y,d,w, type, exclude = NULL)
alpha(obj, mz, mc, k = 98, zcmodel, Tc, Tz)
test.unif.ct(z,y,d,w, type, exclude = 0.05)
findagoodestimator(dat,w)
quickpred(obj, y, type = c("PDF", "SF"))
trans(z,y,d,w,type)

alpha.bisec(theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE)
alpha.bisec.out(A, theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE)
ks(u,v,id,w1,w2, K = 25)

\method{print}{iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{terms}{iqr}(x, \ldots)
\method{model.matrix}{iqr}(object, \ldots)
\method{vcov}{iqr}(object, \ldots)
\method{nobs}{iqr}(object, \ldots)


predict.iqrL.internal(object, level, type = c("coef", "CDF", "QF", "sim"), 
  newdata, p, se = FALSE, \ldots)
\method{print}{iqrL}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.iqrL}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{terms}{iqrL}(x, \ldots)
\method{model.matrix}{iqrL}(object, \ldots)
\method{vcov}{iqrL}(object, \ldots)
\method{nobs}{iqrL}(object, \ldots)

qc.penalty(theta, X, bfun, lambda, pen, H)
fixqc(fit, V, bfun, s, type, tol, maxit, safeit, eps0, 
  lambda, r, maxTry, trace, count, pcross = NULL)
\method{print}{qc.iqr}(x, \ldots)
}
\keyword{internal}

