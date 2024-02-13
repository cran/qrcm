\name{internals}
\alias{internals}

\alias{check.in.iqr}
\alias{check.in0.iqrL}
\alias{check.in.iqrL}
\alias{check.out}
\alias{start.iqr}
\alias{start.iqrL}

\alias{ctiqr.internal}
\alias{iqrL.internal}
\alias{iqrL.fit}

\alias{iqr.ee}
\alias{ciqr.ee}
\alias{ctiqr.ee}
\alias{iciqr.ee}
  
\alias{iobjfun}
\alias{iobjfun.ct}
\alias{iobjfun.ic}

\alias{cov.fun.iqr}

\alias{iqrL.ee}
\alias{cov.fun.iqrL}

\alias{iqr.newton}
\alias{divide.et.impera}
\alias{iqrL.newton}

\alias{pmax0}
\alias{maxind}
\alias{num.fun}
\alias{make.bfun}
\alias{apply_bfun}
\alias{p.bisec}
\alias{p.bisec.internal}
\alias{slp.basis}
\alias{is.slp}
\alias{safesolve}
\alias{middlepoint}

\alias{iqr.waldtest}
\alias{extract.p}
\alias{pred.beta}


\alias{km}
\alias{alpha}
\alias{test.unif.ct}
\alias{findagoodestimator}
\alias{quickpred}
\alias{trans}
\alias{fitgamma}

\alias{alpha.bisec}
\alias{alpha.bisec.out}
\alias{ks}

\alias{print.iqr}
\alias{print.summary.iqr}
\alias{terms.iqr}
\alias{model.matrix.iqr}
\alias{vcov.iqr}
\alias{nobs.iqr}


\alias{predict_iqrL.internal}
\alias{print.iqrL}
\alias{print.summary.iqrL}
\alias{terms.iqrL}
\alias{model.matrix.iqrL}
\alias{vcov.iqrL}
\alias{nobs.iqrL}

\alias{qc.penalty}
\alias{fixqc}
\alias{print.qc.iqr}

\title{Internal Functions}
\description{
Functions for internal use only, or not yet documented.
}
\usage{

check.in.iqr(mf, formula.p, s)
check.in0.iqrL(mf1, mf2)
check.in.iqrL(mf, y,w,formula.p, s)
check.out(theta, S, covar)
start.iqr(y,z,d, x, weights, bfun, df, yy, zz, s, type = "ctiqr")
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
iciqr.ee(theta, y, z, d, X, Xw, bfun, p.star.y, p.star.z, 
  J = TRUE, G, i = FALSE, lambda = 0)
  
iobjfun(theta, y, X, weights, bfun, p.star)
iobjfun.ct(theta, z,y,d,X,weights, bfun, py, pz, type)
iobjfun.ic(fit, V, bfun)

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
middlepoint(y)

iqr.waldtest(obj)
extract.p(model, p, cov = FALSE)
pred.beta(model, p, se = FALSE)


km(z,y,d,w, type, exclude = NULL)
alpha(obj, mz, mc, k = 98, zcmodel, Tc, Tz)
test.unif.ct(z,y,d,w, type, exclude = 0.05)
findagoodestimator(dat, w, type = "ctiqr")
quickpred(obj, y, type = c("PDF", "SF"))
trans(z,y,d,w,type)
fitgamma(y,X,w)

alpha.bisec(theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE)
alpha.bisec.out(A, theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE)
ks(u,v,id,w1,w2, K = 25)

\method{print}{iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{print}{summary.iqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{terms}{iqr}(x, \ldots)
\method{model.matrix}{iqr}(object, \ldots)
\method{vcov}{iqr}(object, \ldots)
\method{nobs}{iqr}(object, \ldots)


predict_iqrL.internal(object, level, type = c("coef", "CDF", "QF", "sim"), 
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

