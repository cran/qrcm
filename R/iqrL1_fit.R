

#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# iqrL functions #############################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
iqrL <- function(fx, fu = ~ slp(u,3), fz = ~ 1, fv = ~ -1 + I(qnorm(v)), id, weights, s.theta, s.phi, data, tol = 1e-5, maxit){

	getmf <- function(cl, mf, f){
		cl$formula <- mf$formula <- f
		m <- match(c("formula", "data", "weights", "id"), names(mf), 0L)
		mf <- mf[c(1L, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1L]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		mf
	}

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	mf1 <- getmf(cl,mf,fx)
	mf2 <- getmf(cl,mf,fz)

	o1 <- attr(mf1, "na.action")
	o2 <- attr(mf2, "na.action")
	o <- unique(c(o1,o2))
	if(!is.null(o)){
		r1 <- attr(mf1, "row.names")
		r2 <- attr(mf2, "row.names")	
		if(!is.null(o1)){mf2 <- mf2[!(r2 %in% names(o1)),, drop = FALSE]}
		if(!is.null(o2)){mf1 <- mf1[!(r1 %in% names(o2)),, drop = FALSE]}
	}

	iqrL.internal(mf1 = mf1, mf2 = mf2, cl = cl, 
		fu = fu, fv = fv,
		s.theta = s.theta, s.phi = s.phi,
		tol = tol, maxit = maxit)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# handles weights, id, and y
# The weight of each individual or cluster is the average of the id's weights.
check.in0.iqrL <- function(mf1, mf2){

	# weights

	if(any((w1 <- model.weights(mf1)) < 0)){stop("negative 'weights'")}
	if(is.null(w1)){w1 <- rep.int(1, nrow(mf1)); alarm <- FALSE}
	else{
	  alarm <- (w1 == 0)
	  sel <- which(!alarm)
	  mf1 <- mf1[sel,, drop = FALSE]
	  mf2 <- mf2[sel,, drop = FALSE]
	  w1 <- w1[sel]
	  w1 <- w1/mean(w1)
	}
	if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
	if((n <- nrow(mf1)) == 0){stop("zero non-NA cases", call. = FALSE)}
	
	# id
	
	if(!any(names(mf1) == "(id)")){stop("'id' is missing")}
	id0 <- mf1[,"(id)"]; id <- as.numeric(as.factor(id0))
	first.id <- which(!duplicated(id))
	u.id0 <- id0[first.id]; u.id <- id[first.id]
	if((n.id <- length(u.id)) == n){stop("all id values are unique", call. = FALSE)}

	# mf2

	mf2.id <- which(names(mf2) == "(id)")
	mf2.w <- which(names(mf2) == "(weights)")
	mf2.z <- (1:ncol(mf2))[-c(mf2.id,mf2.w)]
	if(any(mf2.z)){
		test.mf2 <- mf2[,mf2.id][!duplicated(mf2[,c(mf2.z,mf2.id)])]
		if(length(test.mf2) != n.id){stop("Some of the variables in 'fz' is not at cluster level")}
	}
	mf2 <- mf2[first.id,, drop = FALSE]

	# sorting mf1 and mf2 by id

	o1 <- order(id, mf1[,1])
	mf1 <- mf1[o1,, drop = FALSE]
	id <- id[o1]
	id0 <- id0[o1]
	w1 <- w1[o1]
	
	o2 <- order(u.id)
	mf2 <- mf2[o2,, drop = FALSE]
	u.id <- u.id[o2]
	u.id0 <- u.id0[o2]
	rownames(mf2) <- 1:nrow(mf2)
	
	# output

	list(mf1 = mf1, mf2 = mf2, 
		id0 = id0, id = id, u.id = u.id, u.id0 = u.id0,
		w1 = w1, w2 = c(tapply(w1,id,mean)),
		y = model.response(mf1))
}



#############################################################################################################
#############################################################################################################
#############################################################################################################



check.in.iqrL <- function(mf, y,w,formula.p, s){

	if(!missing(s) && all(s == 0)){stop("'s' cannot be all zero")}
	explore.s <- function(s, dim){
		if(dim == 2){s <- t(s)}
		out <- 1
		if((r <- nrow(s)) > 1){
			for(j in 2:r){
				done <- FALSE; rj <- s[j,]
				for(h in 1:(j - 1)){
					if(all(rj == s[h,])){out[j] <- out[h]; done <- TRUE}
				}
				if(!done){out[j] <- max(out) + 1}
			}
		}
		out
	}
	
	# y, x and b(p)
	
	n <- nrow(mf)
	X <- model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
	termlabelsX <- attr(attr(mf, "terms"), "term.labels")
	assignX <- attr(X, "assign")
	coefnamesX <- colnames(X)
	
	# p1 is used to evaluate the splinefuns. A non-evenly spaced grid, with more values on the tails.
	# p2 is for external use (p.bisec). A grid with p reachable by bisection on the p scale.
	# p3 is for internal use (p.bisec.internal). A grid with p reachable by bisection on the index scale.
	# p4 is like p3, but longer.
	p1 <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = 1000),2,2)
	p2 <- (1:1023)/1024
	p3 <- pbeta(seq.int(qbeta(1/(1000*n),2.5,2.5), qbeta(1 - 1/(1000*n),2.5,2.5), length.out = 1023),2.5,2.5)
	p4 <- pbeta(seq.int(qbeta(1/(1000*n),2.5,2.5), qbeta(1 - 1/(1000*n),2.5,2.5), length.out = 4095),2.5,2.5)
	
	if((use.slp <- is.slp(formula.p))){
		k <- attr(use.slp, "k")
		intercept <- attr(use.slp, "intercept") 	# slp(0) = 0?
		intB <- attr(use.slp, "intB")			# b(p) includes 1?
		assignB <- (1 - intB):k
		termlabelsB <- paste("slp", 1:k, sep = "")
		coefnamesB <- (if(intB) c("(Intercept)", termlabelsB) else termlabelsB)
		k <- k + intB
	}
	else{
		tau <- c(p1,p2,p3,p4)
		B <- model.matrix(formula.p, data = data.frame(p = tau, u = tau, v = tau))
		B1 <- B[1:1000,, drop = FALSE]
		B2 <- B[1001:2023,, drop = FALSE]
		B3 <- B[2024:3046,, drop = FALSE]
		B4 <- B[3047:7141,, drop = FALSE]

		k <- ncol(B)
		assignB <- attr(B, "assign")
		termlabelsB <- attr(terms(formula.p), "term.labels")
		coefnamesB <- colnames(B)
	}
	if(missing(s)){s <- matrix(1,q,k)}
	else{
		if(any(dim(s) != c(q,k))){stop("wrong size of 's.theta' or 's.phi'")}
		if(any(s != 0 & s != 1)){stop("'s.theta' and 's.phi' can only contain 0 and 1")}
	}

	# x singularities (set s = 0 where singularities occur)
	# x is dropped as in a linear model, irrespective of s.

	vx <- qr(X); selx <- vx$pivot[1:vx$rank]
	if(vx$rank < q){s[-selx,] <- 0}

	# b(p) singularities. Dropped row by row, based on s

	if(!use.slp && qr(B3)$rank < k){
		u <- explore.s(s,1)
		for(j in unique(u)){
			sel <- which(s[which(u == j)[1],] == 1)
			if(length(sel) > 1){
				vbj <- qr(B3[,sel, drop = FALSE])
				if((rj <- vbj$rank) < length(sel)){
					s[u == j, sel[-vbj$pivot[1:rj]]] <- 0
				}
			}
		}
	}
	
	# location-scale statistics for x, b(p), and y

	ry <- range(y); my <- ry[1]; My <- ry[2]

	sX <- apply(X,2,sd); mX <- colMeans(X)
	intX <- (length((constX <- which(sX == 0 & mX != 0))) > 0)
	varsX <- which(sX > 0); zeroX <- which(sX == 0 & mX == 0)
	sX[constX] <- X[1,constX]; mX[constX] <- 0; sX[zeroX] <- 1
	if(length(constX) > 1){zeroX <- c(zeroX, constX[-1]); constX <- constX[1]}

	if(!use.slp){
		sB <- apply(B3,2,sd); mB <- colMeans(B3)
		intB <- (length((constB <- which(sB == 0 & mB != 0))) > 0); varsB <- which(sB > 0)
		if(length(varsB) == 0){stop("'fu' and 'fv' must depend on p")}
		if(length(constB) > 1){stop("remove multiple constant functions from 'f.u' or 'f.v'")}
		if(any(sB == 0 & mB == 0)){stop("remove zero functions from 'f.u' or 'f.v'")}
		sB[constB] <- B3[1,constB]; mB[constB] <- 0
	}
	else{
		sB <- rep.int(1, k); mB <- rep.int(0, k)
		if(intB){constB <- 1; varsB <- 2:k}
		else{constB <- integer(0); varsB <- 1:k}
	}

	if(all(s[, varsB] == 0)){stop("the quantile function must depend on p (wrong specification of 's.theta' or 's.phi')")}
	if(!(theta00 <- ((intX & intB) && s[constX, constB] == 1)))
		{my <- 0; My <- sd(y)*5; mX <- rep.int(0,q)}
	else{for(j in varsX){if(any(s[j,] > s[constX,])){mX[j] <- 0}}}
	if(!intB | (intB && any(s[,constB] == 0))){mB <- rep.int(0,k)}

	# Create bfun (only used by post-estimation functions)

	if(!use.slp){
		bfun <- list()
		if(intB){bfun[[constB]] <- function(p, deriv = 0){rep.int(1 - deriv, length(p))}}
		for(j in varsB){bfun[[j]] <- make.bfun(p1,B1[,j])}
		names(bfun) <- coefnamesB
		attr(bfun, "k") <- k
	}
	else{
		bfun <- slp.basis(k - intB, intercept)
		if(!intB){bfun$a[1,1] <- bfun$A[1,1] <- bfun$AA[1,1] <- 0}
		attr(bfun, "intB") <- intB
		B2 <- apply_bfun(bfun,p2, "bfun")
	}

	attr(bfun, "bp") <- B2
	attr(bfun, "p") <- p2

	# first scaling of x, b(p), y

	U <- list(X = X, y = y)
	X <- scale(X, center = mX, scale = sX)
	y <- (y - my)/(My - my)*10
	if(!use.slp){
	  B3 <- scale(B3, center = mB, scale = sB)
	  B4 <- scale(B4, center = mB, scale = sB)
	 }

	# principal component rotations that I can apply to x and b(p); second scaling

	rotX <- (diag(1,q))
	MX <- rep.int(0,q); SX <- rep.int(1,q)
	if(length(varsX) > 0){
		uX <- explore.s(s,1)
		X_in <- rowSums(s)	
		for(j in unique(uX)){
			sel <- which(uX == j)
			if(intX){sel <- sel[sel != constX]}
			if(length(sel) > 1 && X_in[sel[1]] != 0){
				PC <- prcomp(X[,sel], center = FALSE, scale. = FALSE)
				X[,sel] <- PC$x
				rotX[sel,sel] <- PC$rotation
			}
		}
		MX <- colMeans(X); MX[mX == 0] <- 0
		SX <- apply(X,2,sd); SX[constX] <- 1; SX[zeroX] <- 1
		X <- scale(X, center = MX, scale = SX)
	}

	rotB <- (diag(1,k))
	MB <- rep.int(0,k); SB <- rep.int(1,k)
	if(!use.slp){
		uB <- explore.s(s,2)
		B_in <- colSums(s)	
		for(j in unique(uB)){
			sel <- which(uB == j)
			if(intB){sel <- sel[sel != constB]}
			if(length(sel) > 1 && B_in[sel[1]] != 0){
				PC <- prcomp(B3[,sel], center = FALSE, scale. = FALSE)
				B3[,sel] <- PC$x
				rotB[sel,sel] <- PC$rotation
				B4 <- B4%*%rotB
			}
		}
		MB <- colMeans(B3); MB[mB == 0] <- 0
		SB <- apply(B3,2,sd); SB[constB] <- 1
		B3 <- scale(B3, center = MB, scale = SB)
		B4 <- scale(B4, center = MB, scale = SB)
	}

	# Create a pre-evaluated basis (only used internally)

	p <- p3; p_long <- p4
	
	if(!use.slp){
	  bp <- B3; bp_long <- B4
	  dp <- p[-1] - p[-1023]
	  b1p <- num.fun(dp,bp, "der")
	  Bp <- num.fun(dp,bp, "int")
	  BB1 <- num.fun(dp,Bp, "int")[1023,]
	}
	else{
	  k <- attr(bfun, "k")
	  pp <- matrix(, 1023, k + 1); pp_long <- matrix(, 4095, k + 1)
	  pp[,1] <- 1; pp[,2] <- p; pp_long[,1] <- 1; pp_long[,2] <- p_long
	  if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p; pp_long[,j + 1] <- pp_long[,j]*p_long}}
	  bp <- tcrossprod(pp, t(bfun$a)); bp_long <- tcrossprod(pp_long, t(bfun$a))
	  b1p <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))
	  pp <- cbind(pp, pp[,k + 1]*p, pp[,k + 1]*p^2)
	  Bp <- tcrossprod(pp[,2:(k + 2)], t(bfun$A))
	  BB1 <- colSums(bfun$AA)
	    
	  if(!intB){
	    bp <- bp[,-1, drop = FALSE]; bp_long <- bp_long[,-1, drop = FALSE]
	    b1p <- b1p[,-1, drop = FALSE]
	    Bp <- Bp[,-1, drop = FALSE]
	    BB1 <- BB1[-1]
	  }
	}
	BB1 <- matrix(rep(BB1, each = 1023), 1023)
	bpij <- NULL
	for(i in 1:ncol(bp)) {
		bpij <- cbind(bpij, bp * bp[, i])
	}

	# remark: Bp is already BB1 - Bp
	internal.bfun <- list(p = p, bp = bp, b1p = b1p, Bp = BB1 - Bp, bpij = bpij, p_long = p_long, bp_long = bp_long) 
	attr(internal.bfun, "pfun") <- approxfun(c(p[1], 0.5*(p[-1023] + p[-1])), p, method = "constant", rule = 2)


	# output. U = the original variables. V = the scaled/rotated variables.
	# stats.B, stats.X, stats.y = lists with the values use to scale/rotate

	stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
		intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
	stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
		intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
	stats.y <- list(m = my, M = My)

	V <- list(X = X, Xw = X*w, y = y, w = w)
	list(mf = mf, U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y, 
		internal.bfun = internal.bfun, bfun = bfun, s = s)
}





#############################################################################################################
#############################################################################################################
#############################################################################################################


iqrL.internal <- function(mf1, mf2, cl, fu,fv, s.theta, s.phi, tol = 1e-5, maxit){

	## Check-in ####################################

	A0 <- check.in0.iqrL(mf1,mf2)
	mf1 <- A0$mf1; mf2 <- A0$mf2
	id <- A0$id

	A1 <- check.in.iqrL(mf1, A0$y, A0$w1, fu, s.theta)
	V1 <- A1$V; U1 <- A1$U; s.theta <- A1$s
	S1 <- list(B = A1$stats.B, X = A1$stats.X, y = A1$stats.y)
	attributes(A1$bfun) <- c(attributes(A1$bfun), S1$B)
	bfun1 <- A1$internal.bfun

	A2 <- check.in.iqrL(mf2, A0$y, A0$w2, fv, s.phi)
	V2 <- A2$V; U2 <- A2$U; s.phi <- A2$s
	S2 <- list(B = A2$stats.B, X = A2$stats.X, y = A2$stats.y)
	attributes(A2$bfun) <- c(attributes(A2$bfun), S2$B)
	bfun2 <- A2$internal.bfun

	theta00 <- (S1$X$intercept & S1$B$intercept)
	phi00 <- (S2$X$intercept & S2$B$intercept)
	if((theta00 & phi00) && (s.theta[S1$X$const,S1$B$const] & s.phi[S2$X$const, S2$B$const]))
	  {stop("unidentified model, please remove intercept from either 'fu' or 'fv'")}
	if(!(theta00 | phi00)){y <- V1$y}
	else if(theta00){y <- V1$y; S2$y$m <- 0; S2$y$M <- S1$y$M - S1$y$m}
	else{y <- V2$y; S1$y$m <- 0; S1$y$M <- S2$y$M - S2$y$m}

	init <- start.iqrL(y,V1$X,V2$X,id, A0$w1, A0$w2, bfun1,bfun2, s.theta,s.phi, S1,S2)
	theta <- init$theta; phi <- init$phi; alpha <- init$alpha

	r1 <- sum(s.theta)
	r2 <- sum(s.phi)
	npar <- r1 + r2
	n <- nrow(mf1); n.id <- nrow(mf2)

	#################################################

	Fit <- NULL
	fit.ok <- FALSE

	if(missing(maxit)){maxit <- 50 + 5*npar} # maximum n. of iterations
	maxit.theta <- 5 + 2*r1; maxit.phi <- 5 + 2*r2
	safeit.theta <- 2 + r1; safeit.phi <- 2 + r2
	eps <- 0.1 # tuning parameter for optimization
	gTol <- 1
	try.count <- 0
	while(!fit.ok){

		try.count <- try.count + 1

		fit <- iqrL.fit(theta,phi, y,alpha, V1$X,V1$Xw,V2$X,V2$Xw, id, A0$w1,A0$w2, bfun1,bfun2, s.theta,s.phi, 
			maxit.theta, safeit.theta, maxit.phi, safeit.phi, eps, tol, maxit)

		fit.ok <- (fit$fullrank & max(abs(fit$g)/sqrt(n)) < gTol)

		if(fit.ok){
			covar <- try(cov.fun.iqrL(fit, V1$X,V1$Xw,V2$X,V2$Xw, id, A0$w1,A0$w2, 
				bfun1,bfun2, s.theta, s.phi), silent = FALSE)
			fit.ok <- !inherits(covar, "try-error")
		}

		covar.ok <- (if(fit.ok) TRUE else FALSE)
		  # The above is only for consistency with iqr. 
		  # In iqrL, not invertible covar are very frequent, 
		  # and I "force" chol(covar) to be ok using safesolve
		if(fit.ok & covar.ok){break}
		else if(fit.ok){Fit <- fit; Cov <- covar}
		
		if(try.count > 10 && !is.null(Fit)){break}
		if(try.count == 20){break}

		maxit <- maxit + 2
		eps <- eps/2
		safeit.theta <- safeit.theta + 1
		maxit.theta <- maxit.theta + 1
		safeit.phi <- safeit.phi + 1
		maxit.phi <- maxit.phi + 1
		gTol <- gTol * 2
	}

	if(!fit.ok && is.null(Fit)){stop("unable to fit the model: this can be due to severe misspecification")}
	if(!covar.ok){warning("the estimated covariance matrix is deemed to be singular")}
	if(!fit.ok){fit <- Fit; covar <- Cov}
	if(!fit$converged){warning("the algorithm did not converge")}

	# id,y, alpha, y - alpha, u,v, loss

	alpha <- fit$A$alpha*(S2$y$M - S2$y$m)/10 + S2$y$m
	y_alpha <- fit$A$y_alpha*(S1$y$M - S1$y$m)/10 + S1$y$m
	v <- fit$A$v; u <- fit$A$u
	fit.alpha <- data.frame(id = A0$id0, y = A0$y, alpha = alpha[id], 
		y_alpha = y_alpha, v = v[id], u = u)
	loss.theta <- sum(fit$fit.theta$loss*A0$w1)*(S1$y$M - S1$y$m)/10
	loss.phi <- sum(fit$fit.phi$loss*A0$w2)*(S2$y$M - S2$y$m)/10
	loss <- c(obj = loss.theta + loss.phi, obj1 = loss.theta, obj2 = loss.phi)
	attr(loss, "df") <- c(df = sum(s.theta) + sum(s.phi), df1 = sum(s.theta), df2 = sum(s.phi))
	
	# output

	attr(mf1, "assign") <- S1$X$assign
	attr(mf1, "stats") <- S1
	attr(mf1, "all.vars") <- V1
	attr(mf1, "all.vars.unscaled") <- U1
	attr(mf1, "internal.bfun") <- bfun1
	attr(mf1, "bfun") <- A1$bfun
	attr(mf1, "theta") <- fit$fit.theta$par

	mf2 <- delete.response(mf2)
	attr(mf2, "assign") <- S2$X$assign
	attr(mf2, "stats") <- S2
	attr(mf2, "all.vars") <- V2
	attr(mf2, "all.vars.unscaled") <- U2
	attr(mf2, "internal.bfun") <- bfun2
	attr(mf2, "bfun") <- A2$bfun
	attr(mf2, "phi") <- fit$fit.phi$par

	out.theta <- check.out(fit$fit.theta$par, S1, covar$Q.theta)
	out.phi <- check.out(fit$fit.phi$par, S2, covar$Q.phi)

	fit <- list(theta = out.theta$theta, phi = out.phi$theta, obj.function = loss,
		call = cl, converged = fit$converged, n.it = fit$n.it, 
		covar.theta = out.theta$covar, covar.phi = out.phi$covar, 
		mf.theta = mf1, mf.phi = mf2, s.theta = s.theta, s.phi = s.phi,
		fit = fit.alpha
	)

	jnames <- function(x,y){paste(x,y, sep = ":")}
	jnames.theta <- c(sapply(attr(A1$bfun, "coef.names"), jnames, y = S1$X$coef.names))
	jnames.phi <- c(sapply(attr(A2$bfun, "coef.names"), jnames, y = S2$X$coef.names))
	dimnames(fit$covar.theta) <- list(jnames.theta, jnames.theta)
	dimnames(fit$covar.phi) <- list(jnames.phi, jnames.phi)
	
	dimnames(fit$theta) <- dimnames(fit$s.theta) <- list(S1$X$coef.names, S1$B$coef.names)
	dimnames(fit$phi) <- dimnames(fit$s.phi) <- list(S2$X$coef.names, S2$B$coef.names)

	# finish

	class(fit) <- "iqrL"
	fit
}



#############################################################################################################
#############################################################################################################
#############################################################################################################

iqrL.ee <- function(par, x,xw, bfun, 
	p, g = TRUE, H = TRUE, i = FALSE){

	k <- ncol(par)
	n <- nrow(x)
	if(g){
		B <- bfun$Bp[p,, drop = FALSE]
		if(!i){g <- c(crossprod(xw,B))}
		else{g <- NULL; for(j in 1:k){g <- cbind(g,x*B[,j])}}
	}
	if(!H){return(g)}
 

	b1 <- bfun$b1p[p,, drop = FALSE]
	bij <- bfun$bpij[p,, drop = FALSE]

	f <- 1/c(.rowSums(tcrossprod(x, t(par))*b1, n,k)) # pdf
	f <- pmax0(f); f[attr(p, "out")] <- 0
	xwf <- xw*f # weighted x, multiplied by the pdf

	H <- NULL
	count <- 0
	for(i1 in 1:k){
		H.temp <- NULL
		for(i2 in 1:k){
			count <- count + 1
			H.temp <- cbind(H.temp, crossprod(xwf, x*bij[,count]))
		}
		H <- rbind(H, H.temp)
	}
	list(g = g, H = H, f = f, xwf = xwf)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################

start.iqrL <- function(y,x,z,id, w1,w2, bfun1,bfun2, s.theta,s.phi, S1,S2){

	alpha <- c(tapply(y,id,median))
	y_alpha <- y - alpha[id]
	n <- length(y_alpha)
	n.id <- length(alpha)
	pch.fit.ct <- getFromNamespace("pch.fit.ct", ns = "pch")
	predF.pch <- getFromNamespace("predF.pch", ns = "pch")

	##############################################################

	u <- (rank(y_alpha) - 0.5)/n
	if((q1 <- length(S1$X$vars)) > 0){
		df <- max(5, min(15, round(n/30/(q1 + 1))))	
		m1 <- suppressWarnings(pch.fit.ct(y = -log(1 - u), 
			x = cbind(1,x), w = w1, breaks = df))
		u <- 1 - predF.pch(m1)[,3]
	}

	##############################################################

	v <- (rank(alpha) - 0.5)/n.id
	if((q2 <- length(S2$X$vars)) > 0){
		df <- max(5, min(15, round(n.id/30/(q2 + 1))))
		m2 <- suppressWarnings(pch.fit.ct(y = -log(1 - v), 
			x = cbind(1,z), w = w2, breaks = df))
		v <- 1 - predF.pch(m2)[,3]
	}

	##############################################################

	pfun1 <- attr(bfun1, "pfun")
	u <- pfun1(u)
	b1 <- bfun1$bp[match(u, bfun1$p),]
	X <- model.matrix(~ -1 + x:b1)
	X <- X[, c(s.theta) == 1, drop = FALSE]

	pfun2 <- attr(bfun2, "pfun")
	v <- pfun2(v)
	b2 <- bfun2$bp[match(v, bfun2$p),]
	Z <- model.matrix(~ -1 + z:b2)
	Z <- Z[, c(s.phi) == 1, drop = FALSE]
	Z <- Z[id,, drop = FALSE]

	start.ok <- FALSE
	while(!start.ok){
		fit0 <- lm.wfit(cbind(X,Z), y, w1)
		if(any(is.na(fit0$coef))){stop("The specified model is not identifiable. Probably one
		                               covariate is used in both level 1 and 2, with beta(u) and gamma(v)
		                               both having an intercept parameter.")}
		res <- fit0$residuals

		start.ok <- all(res.ok <- (abs(res)/sd(res) < 5))
		if(start.ok | sum(res.ok) < 0.5*n){break}
		X <- X[res.ok,, drop = FALSE]
		Z <- Z[res.ok,, drop = FALSE]
		y <- y[res.ok]
		w1 <- w1[res.ok]
	}

	##############################################################

	rr1 <- length(s.theta)
	r1 <- sum(s.theta)
	theta <- rep.int(0, rr1)
	theta[c(s.theta) == 1] <- fit0$coef[1:r1]
	theta <- matrix(theta, ncol(x))

	rr2 <- length(s.phi)
	r2 <- sum(s.phi)
	phi <- rep.int(0, rr2)
	phi[c(s.phi) == 1] <- fit0$coef[(r1 + 1):(r1 + r2)]
	phi <- matrix(phi, ncol(z))

	##############################################################

	list(theta = theta, phi = phi, alpha = .rowSums(tcrossprod(z, t(phi))*b2, nrow(z),ncol(phi)))
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



# Important remark: when \hat alpha_i is on boundary (v_i = 0 or 1)
# and f(alpha) = 0, G(alpha_i) != 0 and standard asymptotics fail.
# In an earlier version of this code, I artificially "removed" these observations 
# (assuming they are not "too many") by letting H(csi, alpha_i) = O(csi, alpha_i) = 0, 
# which imposes cov(\hat alpha_i, \hat csi) = 0. In the final version I don't do it,
# but I left the old commands (with comments).

cov.fun.iqrL <- function(fit, x,xw,z,zw, id, w1,w2, bfun1,bfun2, s.theta, s.phi){

	# notation:

	# H = hessian
	# V = inverse of H
	# G = gradient
	# O = outer product

	# 1 = theta
	# 2 = phi
	# csi = (theta,phi)
	# a = alpha

	# N = tot. n. of obs
	# n = n. of distinct id
	# m = cluster sizes

	################################################################################
	################################################################################

	N <- nrow(x); n <- nrow(z); m <- c(table(id))
	pp <- mean(m <= 3) # proportion of very small clusters
	lambda <- 1 - 0.1*(pp >= 0.1) - 0.1*(pp >= 0.25) # lambda for safesolve

	theta <- fit$fit.theta$par
	phi <- fit$fit.phi$par
	q1 <- nrow(theta); k1 <- ncol(theta)
	q2 <- nrow(phi); k2 <- ncol(phi)
	r1 <- length(i1 <- which(s.theta == 1)); rr1 <- q1*k1
	r2 <- length(i2 <- which(s.phi == 1)); rr2 <- q2*k2
	ind1 <- 1:r1; ind2 <- (r1 + 1):(r1 + r2)

	p1 <- fit$A$p1
	p2 <- fit$A$p2
	u <- fit$A$u
	v <- fit$A$v
	f1 <- fit$fit.theta$f#; f1pos <- (f1 > 0) # f(y - alpha)
	f2 <- fit$fit.phi$f#; f2pos <- (f2 > 0) # f(alpha)
	xwf1 <- fit$fit.theta$xwf # w1*x*f(y - alpha)
	zwf2 <- fit$fit.phi$xwf # w2*z*f(alpha)

	################################################################################
	################################################################################
	############################## Hessian #########################################
	################################################################################
	################################################################################

	# H(csi,csi) ###################################################################
	
	H.csi <- matrix(0, r1 + r2, r1 + r2)
	H.csi[ind1,ind1] <- fit$fit.theta$H/N # theta
	H.csi[ind2,ind2] <- fit$fit.phi$H/n # phi

	# H(alpha,theta) and H(alpha,phi) ##############################################

	bu <- bfun1$bp[p1,, drop = FALSE]
	bv <- bfun2$bp[p2,, drop = FALSE]
	H.a.1 <- matrix(,n,r1)
	H.a.2 <- matrix(,n,r2)
	for(i in 1:n){
		ind <- which(id == i)
		H.a.1[i,] <- c(crossprod(xwf1[ind,,drop = FALSE],bu[ind,, drop = FALSE]))[i1]
		H.a.2[i,] <- c(crossprod(t(zwf2[i,]),t(bv[i,])))[i2]
	}
	H.a.csi <- cbind(H.a.1/sqrt(N*m),-H.a.2/sqrt(n*m))
	#H.a.csi <- H.a.csi*f2pos # artificially remove where \hat(alpha) is on boundary
	# NOTE: in H.a.csi, values that may cause problems are NOT NECESSARILY the very large ones.

	# H(alpha,alpha) ###############################################################

	H.a.a <- c(tapply(w1*f1,id,sum))/m + w2*f2/m # diagonal only
	#H.a.a[!f2pos] <- 1 # an arbitrary value, to avoid V.a.a = Inf
	H.a.a[H.a.a == 0] <- 1e-3 # very rare issue
	
	################################################################################
	################################################################################
	########################## Inverse of H ########################################
	################################################################################
	################################################################################

	# Elements of V = H^(-1), using block inversion #################################
  	# Note: V11 is the only inverse I have to calculate #############################

	V.a.a <- 1/H.a.a
	A <- H.a.csi*V.a.a
	AA <- crossprod(A, H.a.csi)
	V11 <- safesolve(H.csi, AA, lambda); warn1 <- V11$warn; V11 <- V11$inv
	V12 <- -tcrossprod(V11, A); V21 <- t(V12)

	################################################################################
	################################################################################
	########################## Gradient ############################################
	################################################################################
	################################################################################

	G.1 <- iqrL.ee(theta, x,xw, bfun1, p1, g = TRUE, H = FALSE, i = TRUE)[,i1,drop = FALSE]
	G.2 <- iqrL.ee(phi, z,zw, bfun2, p2, g = TRUE, H = FALSE, i = TRUE)[,i2,drop = FALSE]
	G.2.long <- G.2[id,, drop = FALSE]
	G.a.1 <- (0.5 - u)
	G.a.2 <- (v - 0.5)

	################################################################################
	################################################################################
	########################## Outer product #######################################
	################################################################################
	################################################################################

	# O(csi,csi) ###################################################################

	O.csi <- matrix(0, r1 + r2, r1 + r2)
	O.csi[ind1,ind1] <- crossprod(w1*G.1, G.1)/N
	O.csi[ind2,ind2] <- crossprod(w2*G.2, G.2)/n
	
	U12 <- crossprod(w1*G.1, G.2.long)/sqrt(N*n) # note: asymptotically zero
	# I should do this: O.csi[ind1,ind2] <- U12; O.csi[ind2,ind1] <- t(U12) # or just do nothing, because U12 --> 0.
	# I use U12, but this may cause O.csi to be not positive definite. To overcome this, I do:
	O.csi.bis <- matrix(0, r1 + r2, r1 + r2); O.csi.bis[ind1,ind2] <- U12; O.csi.bis[ind2,ind1] <- t(U12)
	# Now O.csi should be equal to O.csi + O.csi.bis, and I compute it with safesolve:
	O.csi <- safesolve(O.csi, -O.csi.bis, lambda)$X

	# O(alpha,csi) #################################################################

	G.a.1w <- w1*G.a.1
	O.csi.a <- matrix(, r1 + r2, n)
	for(j in ind1){
	  O.csi.a[j,] <- c(tapply(G.1[,j]*G.a.1w,id,sum))/sqrt(N*m) + c(tapply(w1*G.1[,j]*G.a.2[id], id, sum))/sqrt(N*m)
	 }
	for(j in ind2){
	  O.csi.a[j,] <- c(tapply(G.2.long[,j - r1]*G.a.1w,id,sum))/sqrt(n*m) + G.2[,j - r1]*G.a.2*w2/sqrt(n*m)
	}
	# O.csi.a[,!f2pos] <- 0 # artificially remove where \hat(alpha) is on boundary

	# O(alpha,alpha) ###############################################################

	O.a.a <- c(tapply(w1*G.a.1^2,id,sum))/m + w2*G.a.2^2/m + 2*c(tapply(w1*G.a.1*G.a.2[id], id, sum))/m # diagonal only
	O.a.a <- pmax0(O.a.a)

	################################################################################
	################################################################################
	########################## Assembling ##########################################
	################################################################################
	################################################################################

	Q1 <- V11%*%O.csi%*%V11
	Q2 <- V12%*%t(O.csi.a)%*%V11
	Q3 <- t(Q2)
	Q4 <- V12%*%(V21*O.a.a)

	Q14 <- Q1 + Q4
	Q23 <- Q2 + Q3

	Q <- safesolve(Q14, -Q23, lambda)
	warn2 <- Q$warn; Q <- Q$X
	
	Q <- (Q + t(Q))/2 # to avoid numerical asymmetry
	Q[ind1,ind1] <- Q[ind1,ind1]/N
	Q[ind2,ind2] <- Q[ind2,ind2]/n
	Q[ind1,ind2] <- Q[ind1,ind2]/sqrt(N*n)
	Q[ind2,ind1] <- Q[ind2,ind1]/sqrt(N*n)
	
	# output

	Q.theta <- matrix(0,rr1,rr1)
	Q.theta[i1,i1] <- Q[ind1,ind1]
	Q.phi <- matrix(0,rr2,rr2)
	Q.phi[i2,i2] <- Q[ind2,ind2]

	# messages
	
	if(warn1 | warn2){warning(
	  "The asymptotic covariance matrix was not well-defined. 
	  Standard errors may be incorrect.", call. = FALSE)}

	list(Q = Q, Q.theta = Q.theta, Q.phi = Q.phi)
}




#############################################################################################################
#############################################################################################################
#############################################################################################################

iqrL.newton <- function(par, y,x,xw, bfun, s, tol, maxit, safeit, eps){
 
  q <- nrow(par)
  k <- ncol(par)
  s <- c(s == 1)
  
  p <- p.bisec.internal(par, y,x, bfun$bp)
  g <- iqrL.ee(par, x,xw, bfun, p, g = TRUE, H = FALSE)[s]
  conv <- FALSE

  # Preliminary safe iterations, only g is used

  for(i in 1:safeit){

    if(conv | max(abs(g)) < tol){break}
    u <- rep.int(0, q*k)
    u[s] <- g
    delta <- matrix(u, q,k)	
    delta[is.na(delta)] <- 0
    cond <- FALSE

    while(!cond){

      new.par <- par - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p <- p.bisec.internal(new.par, y,x, bfun$bp)
      g1 <- iqrL.ee(new.par, x,xw, bfun, p, g = TRUE, H = FALSE)[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }

    if(conv){break}
    g <- g1
    par <- new.par
    eps <- min(eps*2,0.1)
  }
  
  # Newton-Raphson
  
  alg <- "nr"
  conv <- FALSE
  eps <- 0.1
  A <- iqrL.ee(par, x,xw, bfun, p, g = FALSE, H = TRUE)
  H <- A$H[s,s, drop = FALSE]
  H <- H + diag(0.0001, nrow(H))

  for(i in 1:maxit){

    if(conv | max(abs(g)) < tol){break}
    
    ####
    
    H1 <- try(chol(H), silent = TRUE)
    err <- inherits(H1, "try-error")
    if(!err){
      if(alg == "gs"){alg <- "nr"; eps <- 1}
      delta <- chol2inv(H1)%*%g
    }
    else{
      if(alg == "nr"){alg <- "gs"; eps <- 1}
      delta <- g
    }

    u <- rep.int(0, q*k)
    u[s] <- delta
    delta <- matrix(u, q,k)	
    delta[is.na(delta)] <- 0
    cond <- FALSE
    while(!cond){
      new.par <- par - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p <- p.bisec.internal(new.par, y,x, bfun$bp)
      g1 <- iqrL.ee(new.par, x,xw, bfun, p, g = TRUE, H = FALSE)[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }
    
    if(conv){break}
    g <- g1
    par <- new.par
    A <- iqrL.ee(par, x,xw, bfun, p, g = FALSE, H = TRUE)
    H <- A$H[s,s, drop = FALSE]
    if(i > 1){eps <- min(eps*10,1)}
    else{eps <- min(eps*10,0.1)}
  }
  
  p <- p.bisec.internal(par, y,x, bfun$bp)
  Fy <- bfun$p[p]
  B <- bfun$Bp[p,, drop = FALSE]
  loss <- y*(Fy - 0.5) + .rowSums(tcrossprod(x, t(par))*B, length(y),k)

  list(par = par,
      converged = (i < maxit), n.it = i, 
      p = p, F = Fy, f = A$f, xwf = A$xwf,
      g = g, H = H, fullrank = (alg == "nr"),
	loss = loss)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# if long = TRUE, a grid of 4095 elements is used instead of the "standard" grid of 1023.
alpha.bisec <- function(theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE){
  
  n <- nrow(x)
  n.id <- nrow(z)
  k1 <- ncol(theta)
  k2 <- ncol(phi)
  xtheta <- tcrossprod(x, t(theta))
  zphi <- tcrossprod(z, t(phi))
  
  r <- 10 + 2*long; cc <- 2^(r - 1)
  if(long){
    bfun1$bp <- bfun1$bp_long; bfun2$bp <- bfun2$bp_long
    bfun1$p <- bfun1$p_long; bfun2$p <- bfun2$p_long  
  }
  
  v <- as.integer(rep.int(cc,n.id))
  for(i2 in 2:(r + 1)){
    
    alpha <- .rowSums(zphi*bfun2$bp[v,, drop = FALSE], n.id, k2)
    y_alpha <- y - alpha[id]
    u <- as.integer(rep.int(cc,n))
    for(i1 in 2:r){
      eta <- .rowSums(xtheta*bfun1$bp[u,, drop = FALSE], n, k1)
      u <- u + as.integer(sign(y_alpha - eta)*(2^(r - i1)))
    }		
    obj <- tapply(w1*(0.5 - bfun1$p[u]),id,sum) + w2*(bfun2$p[v] - 0.5)
    if(i2 < (r + 1)){v <- v - as.integer(sign(obj)*(2^(r - i2)))}
  }
  list(alpha = alpha, y_alpha = y_alpha, v = v, u = u, obj = obj, xtheta = xtheta, zphi = zphi)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


# Continues to bisect beyond v = 0 and v = 1
alpha.bisec.out <- function(A, theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = FALSE){
  
  r <- 10 + 2*long; cc <- 2^(r - 1)
  outL <- which(A$v == 1)
  outR <- which(A$v == 2^r - 1)
  out <- c(outL, outR)
  if(!any(out)){return(A)}
 
  ###
  
  out_left <- c(rep(TRUE, length(outL)), rep(FALSE, length(outR)))
  out_right <- !out_left
  w <- which(id %in% out)
  n <- length(w); k1 <- ncol(theta)
  xtheta <- A$xtheta[w,, drop = FALSE]
  id <- id[w]; id <- rep.int(1:length(out), table(id))
  y <- y[w]; u0 <- A$u[w]
  w1 <- w1[w]; w2 <- w2[out]
  if(long){bfun1$bp <- bfun1$bp_long; bfun1$p <- bfun1$p_long}

  # internal bound for alpha
  alpha1 <- A$alpha[out]
  
  # external bound for alpha
  u <- as.integer(rep.int(cc,n))
  y_alpha <- .rowSums(xtheta*bfun1$bp[u,, drop = FALSE], n, k1)
  alpha2 <- t(sapply(tapply(y - y_alpha, id, range), I))
  alpha2 <- alpha2[,1]*out_left + alpha2[,2]*out_right
  
  # bisection

  v <- 0.5*(out_right - out_left) # this is actually v - 0.5
  for(j in 1:100){

    alpha <- 0.5*(alpha1 + alpha2)
    y_alpha <- y - alpha[id]
    u <- as.integer(rep.int(cc,n))
    for(i1 in 2:r){
      eta <- .rowSums(xtheta*bfun1$bp[u,, drop = FALSE], n, k1)
      u <- u + as.integer(sign(y_alpha - eta)*(2^(r - i1)))
    }
    obj <- tapply(w1*(0.5 - bfun1$p[u]),id,sum) + w2*v

    if(max(abs(u - u0)) < 2){break}
    u0 <- u
    obj.temp <- obj*out_left - obj*out_right
    opos <- which(obj.temp > 0); alpha1[opos] <- alpha[opos]
    oneg <- which(obj.temp < 0); alpha2[oneg] <- alpha[oneg]
  }
  A$u[w] <- u
  A$obj[out] <- obj
  A$alpha[out] <- alpha
  A$y_alpha[w] <- y_alpha
  A
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


iqrL.fit <- function(theta,phi, y,alpha, x,xw,z,zw, id, w1,w2, bfun1,bfun2, s.theta, s.phi, 
	maxit.theta, safeit.theta, maxit.phi, safeit.phi, eps, tol, maxit){

	par0 <- c(theta,phi); L0 <- Inf; long <- FALSE

	for(i in 1:maxit){

	  A.new <- alpha.bisec(theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = long)
	  A.new <- alpha.bisec.out(A.new, theta,phi,y,x,z,id,w1,w2,bfun1,bfun2, long = long)
	  alpha <- A.new$alpha
    alpha <- round(alpha, 8) # this seems to avoid numerical issues.
    
		fit.theta.new <- iqrL.newton(theta, y - alpha[id],x,xw, bfun1, s.theta, 
			tol = 1e-6, maxit.theta, safeit.theta, eps)

		fit.phi.new <- iqrL.newton(phi, alpha,z,zw, bfun2, s.phi, 
			tol = 1e-6, maxit.phi, safeit.phi, eps)
		
		new.theta <- fit.theta.new$par; new.phi <- fit.phi.new$par
		par <- c(new.theta, new.phi)
		L <- sum(fit.theta.new$loss) + sum(fit.phi.new$loss)

		if(L > L0 & !long){long <- TRUE; next}
		if(i > 1){if(max(abs(par - par0)) < tol | L0 - L < tol){break}}
		A <- A.new; fit.theta <- fit.theta.new; fit.phi <- fit.phi.new
		theta <- new.theta; phi <- new.phi
		par0 <- par; L0 <- L
		
		# Comment: I could update eps: restarting with the last used eps, MULTIPLIED BY SOMETHING
		# (e.g., eps*8, eps*16: otherwise it will stop immediately...) will usually reduce the number
		# of iterations to convergence.
	}

	A$p1 <- fit.theta$p; A$u <- fit.theta$F
	A$p2 <- fit.phi$p; A$v <- fit.phi$F
	list(fit.theta = fit.theta, fit.phi = fit.phi, A = A, 
		n.it = i, converged = (i < maxit), 
		fullrank = (fit.theta$fullrank & fit.phi$fullrank),
		g = c(fit.theta$g, fit.phi$g)
	)
}






