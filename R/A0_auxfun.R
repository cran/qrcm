#' @importFrom stats model.response model.weights model.matrix terms model.frame delete.response coef 
#' @importFrom stats printCoefmat .getXlevels vcov nobs predict formula update
#' @importFrom stats pbeta qbeta dgamma rgamma rexp pchisq qchisq runif pnorm  
#' @importFrom stats nlm prcomp approxfun integrate splinefun sd median lm.wfit weighted.mean 
#' @importFrom survival Surv survfit coxph
#' @importFrom graphics plot points abline polygon
#' @importFrom grDevices adjustcolor
#' @importFrom utils menu setTxtProgressBar txtProgressBar tail getFromNamespace
#' @importFrom icenReg ic_np getFitEsts
#' @import pch




# VARIOUS AUXILIARY FUNCTIONS, which are used by both iqr and iqrL.

#############################################################################################################
#############################################################################################################
#############################################################################################################


pmax0 <- function(x){(x + abs(x))/2}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# Compute either bfun or b1fun. Note that I never need them both in predictions!
apply_bfun <- function(bfun, p, fun = c("bfun", "b1fun")){
  
  k <- attr(bfun, "k")
  n <- length(p)
  
  if(inherits(bfun, "slp.basis")){
    pp <- matrix(, n, k + 1)
    pp[,1] <- 1; pp[,2] <- p
    if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p}}
    if(fun == "bfun"){out <- tcrossprod(pp, t(bfun$a))}
    else{out <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))}
    if(!attr(bfun, "intB")){out <- out[,-1, drop = FALSE]}
  }
  else{
    out <- matrix(,n,k)
    if(fun == "bfun"){for(j in 1:k){out[,j] <- bfun[[j]](p)}}
    else{for(j in 1:k){out[,j] <- bfun[[j]](p, deriv = 1)}}
  }
  out
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# Bisection for external use. Precision about 1e-6
p.bisec <- function(theta, y,X, bfun, n.it = 20){

  n <- length(y); k <- ncol(theta)
  bp <- attr(bfun, "bp")
  p <- attr(bfun, "p")
  Xtheta <- tcrossprod(X, t(theta))

  eta <- tcrossprod(Xtheta,bp[512,, drop = FALSE])
  m <- as.integer(512 + sign(y - eta)*256)
 
  for(i in 3:10){
    eta <- .rowSums(Xtheta*bp[m,, drop = FALSE], n, k)
    m <- m + as.integer(sign(y - eta)*(2^(10 - i)))
  }
  m <- p[m]
  
  for(i in 11:n.it){
    bp <- apply_bfun(bfun, m, "bfun")
    delta.m <- y - .rowSums(Xtheta*bp, n,k)
    m <- m + sign(delta.m)/2^i
  }
  
  m <- c(m)
  
  out.l <- which(m == 1/2^n.it)
  out.r <- which(m == 1 - 1/2^n.it)
  m[out.l] <- 0; m[out.r] <- 1
  attr(m, "out") <- c(out.l, out.r)
  attr(m, "out.r") <- out.r
  
  m
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# Bisection for internal use. Precision about 0.001
p.bisec.internal <- function(theta, y,X,bp){

	n <- length(y); k <- ncol(theta)

	Xtheta <- tcrossprod(X, t(theta))
	eta <- tcrossprod(Xtheta,bp[512,, drop = FALSE])
	m <- as.integer(512 + sign(y - eta)*256)

	for(i in 3:10){
		eta <- .rowSums(Xtheta*bp[m,, drop = FALSE], n, k)
		m <- m + as.integer(sign(y - eta)*(2^(10 - i)))
	}

	out.l <- which(m == 1)
	out.r <- which(m == 1023)

	attr(m, "out") <- c(out.l, out.r)
	attr(m, "out.r") <- out.r

	m
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



# for external use, only returns b(p)
make.bfun <- function(p,x){
  n <- length(x)
  x1 <- x[1:(n-1)]
  x2 <- x[2:n]
  if(all(x1 < x2) | all(x1 > x2)){method <- "hyman"}
  else{method <- "fmm"}
  splinefun(p,x, method = method)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


num.fun <- function(dx,fx, op = c("int", "der")){
  n <- length(dx) + 1
  k <- ncol(fx)
  fL <- fx[1:(n-1),, drop = FALSE]
  fR <- fx[2:n,, drop = FALSE]
  
  if(op == "int"){out <- apply(rbind(0, 0.5*dx*(fL + fR)),2,cumsum)}
  else{
    out <- (fR - fL)/dx
    out <- rbind(out[1,],out)
  }
  out
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


extract.p <- function(model,p, cov = FALSE){

	theta <- model$coefficients
	v <- model$covar
	q <- nrow(theta)
	k <- ncol(theta)

	bfun <- attr(model$mf, "bfun")
	pred.p <- apply_bfun(bfun, p, "bfun")
	beta <- c(pred.p%*%t(theta))

	cov.beta <- matrix(NA,q,q)
	for(j1 in 1:q){
		w1 <- seq.int(j1,q*k,q)
		for(j2 in j1:q){
			w2 <- seq.int(j2,q*k,q)
			cc <- v[w1,w2, drop = FALSE]
			cov.beta[j1,j2] <- cov.beta[j2,j1] <- pred.p%*%cc%*%t(pred.p)
		}
	}
	se <- sqrt(diag(cov.beta))
	z <- beta/se
	out <- cbind(beta, se, z, 2*pnorm(-abs(z)))
	colnames(out) <- c("Estimate", "std.err", "z value", "p(>|z|))")
	rownames(out) <- colnames(cov.beta) <- rownames(cov.beta) <- rownames(theta)
	if(cov){list(coef = out, cov = cov.beta)}
	else{list(coef = out)}
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


pred.beta <- function(model, p, se = FALSE){

	if(se){
		Beta <- NULL
		SE <- NULL
		for(j in p){
			b <- extract.p(model,j)$coef
			Beta <- rbind(Beta, b[,1])
			SE <- rbind(SE, b[,2])
		}
		out <- list()
		for(j in 1:ncol(Beta)){
			low <- Beta[,j] - 1.96*SE[,j]
			up <- Beta[,j] + 1.96*SE[,j]
			out[[j]] <- data.frame(p = p, beta = Beta[,j], se = SE[,j], low = low, up = up)
		}
		names(out) <- rownames(model$coefficients)
		return(out)
	}
	else{
		theta <- model$coefficients
		beta <- apply_bfun(attr(model$mf, "bfun"), p, "bfun")%*%t(theta)
		out <- list()
		for(j in 1:nrow(theta)){out[[j]] <- data.frame(p = p, beta = beta[,j])}
		names(out) <- rownames(theta)
		return(out)
	}
}




#############################################################################################################
#############################################################################################################
#############################################################################################################


iqr.waldtest <- function(obj){
	bfun <- attr(obj$mf, "bfun")
	ax <- attr(obj$mf, "assign")
	ap <- attr(bfun, "assign")
	theta <- obj$coefficients
	q <- nrow(theta)
	k <- ncol(theta)
	s <- obj$s
	cc <- obj$covar
	ind.x <- rep(ax,k)
	ind.p <- sort(rep.int(ap,q))
	K <- tapply(rowSums(s), ax, sum)
	Q <- tapply(colSums(s), ap, sum)

	testx <- testp <- NULL

	if(q > 1){
		for(i in unique(ax)){
			theta0 <- c(theta[which(ax == i),])
			w <- which(ind.x == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]
			
			if(length(theta0) == 0){tx <- NA}
			else{tx <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testx <- c(testx, tx)
		}
		testx <- cbind(testx, df = K, pchisq(testx, df = K, lower.tail = FALSE))
		colnames(testx) <- c("chi-square", "df", "P(> chi)")

		nx <- attr(attr(obj$mf, "terms"), "term.labels")
		if(attr(attr(obj$mf, "terms"), "intercept") == 1){nx <- c("(Intercept)", nx)}
		rownames(testx) <- nx
	}
	if(k > 1){
		for(i in unique(ap)){
			theta0 <- c(theta[,which(ap == i)])
			w <- which(ind.p == i)
			c0 <- cc[w,w, drop = FALSE]

			theta0 <- theta0[theta0 != 0]
			w <- which(rowSums(c0) != 0)
			c0 <- c0[w,w, drop = FALSE]

			if(length(theta0) == 0){tp <- NA}
			else{tp <- t(theta0)%*%chol2inv(chol(c0))%*%t(t(theta0))}
			testp <- c(testp, tp)
		}
		testp <- cbind(testp, df = Q, pchisq(testp, df = Q, lower.tail = FALSE))
		colnames(testp) <- c("chi-square", "df", "P(> chi)")
		np <- attr(bfun, "term.labels")
		if(any(ap == 0)){np <- c("(Intercept)", np)}
		rownames(testp) <- np
	}

	list(test.x = testx, test.p = testp)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
test.fit <- function (object, ...){UseMethod("test.fit")}


#############################################################################################################
#############################################################################################################
#############################################################################################################

check.out <- function(theta, S, covar){

	blockdiag <- function(A, d, type = 1){
		h <- nrow(A); g <- d/h
		if(type == 1){
			out <- diag(1,d)
			for(j in 1:g){ind <- (j*h - h  + 1):(j*h); out[ind,ind] <- A}
		}
		else{
			out <- matrix(0,d,d)
			for(i1 in 1:h){
				for(i2 in 1:h){
					ind1 <- (i1*g - g  + 1):(i1*g)
					ind2 <- (i2*g - g  + 1):(i2*g)
					out[ind1, ind2] <- diag(A[i1,i2],g)
				}
			}
			out <- t(out)
		}
		out
	}

	mydiag <- function(x){
		if(length(x) > 1){return(diag(x))}
		else{matrix(x,1,1)}
	}

	th <- cbind(c(theta))
	q <- nrow(theta)
	k <- ncol(theta)
	g <- q*k
	aX <- S$X; ay <- S$y; aB <- S$B
	cX <- aX$const; cB <- aB$const

	##########################

	A <- blockdiag(mydiag(1/aX$S), g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$M; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$S),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$M; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(aX$rot,g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	A <- blockdiag(t(aB$rot),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	##########################

	A <- blockdiag(mydiag(1/aX$s),g)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aX$intercept){
		A <- diag(1,q); A[cX,] <- -aX$m/aX$s[cX]; A[cX, cX] <- 1
		A <- blockdiag(A,g)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	A <- blockdiag(mydiag(1/aB$s),g,2)
	th <- A%*%th
	covar <- A%*%covar%*%t(A)

	if(aB$intercept){
		A <- diag(1,k); A[,cB] <- -aB$m/aB$s[cB]; A[cB, cB] <- 1
		A <- blockdiag(A,g,2)
		th <- A%*%th
		covar <- A%*%covar%*%t(A)
	}

	##########################

	v <- (ay$M - ay$m)/10
	th <- th*v
	covar <- covar*(v^2)
	theta <- matrix(th,q,k)
	theta[cX,cB] <- theta[cX,cB] + ay$m/aB$s[cB]/aX$s[cX]
	
	list(theta = theta, covar = covar)
}

#############################################################################################################
#############################################################################################################
#############################################################################################################

# The function below is quite arbitrary, but somehow motivated. And seems to work well.
# Of course, please don't expect meaningful results with asymptotics, if each cluster has size T = 2...

# Goal: to compute the inverse of A - B. Don't use this function for a generic A and B!
 # This function is specially designed for MY estimation problem, in which I know the following facts:
  ### A is always positive definite (H.csi is hopefully ok, O.csi and Q14 are quadratic forms)
  ### B may or may not be positive definite (AA will, O.csi.bis will not, Q23 maybe).

# B incorporates the extra variance associated with estimation of alpha. This means that:
  ### for B = 0, I can solve A - B
  ### B = 0 corresponds to "the individual effects are known" 
  ### asymptotically, B = 0.

# Instead of solving A - B, I solve A - lambda*B, lambda <= 1.
# I have two goals: (1) make sure A - lambda*B can be solved, and (2) controlling the variance.
# This is why, although lambda = 1 may solve the problem, I may start from a smaller lambda.
# This alone may help avoiding very large estimated standard errors. Empirically,
# the standard errors tend to be overestimated, and this is a heuristic fix that,
# in the worst case, does not do too bad.

# The lambda from which to start is defined in cov.fun.iqrL, and works as follows:
  # Define p = % of clusters with <= 3 observations.
  # If p < 0.1, I start from lambda = 1.
  # If 0.1 <= p < 0.25, I start from lambda = 0.9.
  # If p >= 0.25, I start from lambda = 0.8.
  # I don't go further down, because this may cause bias in the estimated standard errors.

safesolve <- function(A,B, lambda){
  
  if(any(eigen(A)$values <= 0)){stop("A is not definite positive")} 
  # this will exit and repeat estimation from scratch
  
  n <- nrow(A)
  lambda <- seq(lambda,0,length = 20)
  
  for(i in 1:length(lambda)){
    X <- A - lambda[i]*B
    inv <- try(chol2inv(chol(X)), silent = TRUE)
    if(!inherits(inv, "try-error")){break}
  }

  list(X = X, inv = inv, lambda = lambda[i], warn = (i > 1))
}

# Indexes of a matix A where A == max(A)
maxind <- function(A){
  
  w <- which.max(A)
  n <- nrow(A)
  m <- ncol(A)
  
  row <- w %% n
  if(row == 0){row <- n}
  col <- floor(w/n) + (row < n)
  
  c(row, col)
}