


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# iqr functions #############################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
iqr <- function(formula, formula.p = ~ slp(p,3), weights, data, s, tol = 1e-6, maxit, remove.qc = FALSE){
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	ctiqr.internal(mf = mf,cl = cl, formula.p = formula.p, tol = tol, maxit = maxit, s = s, remove.qc)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


check.in.iqr <- function(mf, formula.p, s){

	if(!(any(all.vars(formula.p) == "p"))){stop("the quantile function must depend on p")}
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

	# weights

	if(any((weights <- model.weights(mf)) < 0)){stop("negative 'weights'")}
	if(is.null(weights)){weights <- rep.int(1, nrow(mf)); alarm <- FALSE}
	else{
	alarm <- (weights == 0)
	  sel <- which(!alarm)
	  mf <- mf[sel,]
	  weights <- weights[sel]
	  weights <- weights/mean(weights)
	}
	if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
	if((n <- nrow(mf)) == 0){stop("zero non-NA cases", call. = FALSE)}
	
	# y,z,d

	zyd <- model.response(mf)
	type <- attributes(zyd)$type
	zyd <- cbind(zyd)

	if(is.null(type)){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- rep.int(1,n); type <- fittype <- "iqr"}
	else if(type == "right"){
	  y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
	  type <- (if(any(d == 0)) "ciqr" else "iqr")
	  fittype <- "ciqr"
	}
	else if(type == "counting"){
	  z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- fittype <- "ctiqr"
	  if(all(z < min(y))){type <- (if(any(d == 0)) "ciqr" else "iqr")}
	}
	attr(type, "fittype") <- fittype
	if(!(any(d == 1))){stop("all data are censored")}

	# x and b(p)

	X <- model.matrix(attr(mf, "terms"), mf); q <- ncol(X)
	termlabelsX <- attr(attr(mf, "terms"), "term.labels")
	assignX <- attr(X, "assign")
	coefnamesX <- colnames(X)
	
	# p1 is used to evaluate the splinefuns. A non-evenly spaced grid, with more values on the tails.
	# p2 is for external use (p.bisec). A grid with p reachable by bisection on the p scale.
	# p3 is for internal use (p.bisec.internal). A grid with p reachable by bisection on the index scale.
	p1 <- pbeta(seq.int(qbeta(1e-6,2,2), qbeta(1 - 1e-6,2,2), length.out = 1000),2,2)
	p2 <- (1:1023)/1024
	p3 <- pbeta(seq.int(qbeta(1/(1000*n),2.5,2.5), qbeta(1 - 1/(1000*n),2.5,2.5), length.out = 1023),2.5,2.5)

  
	if((use.slp <- is.slp(formula.p))){
		k <- attr(use.slp, "k")
		intercept <- attr(use.slp, "intercept") 	# slp(0) = 0?
		intB <- attr(use.slp, "intB")			        # b(p) includes 1?
		assignB <- (1 - intB):k
		termlabelsB <- paste("slp", 1:k, sep = "")
		coefnamesB <- (if(intB) c("(Intercept)", termlabelsB) else termlabelsB)
		k <- k + intB
	}
	else{
	  B <- model.matrix(formula.p, data = data.frame(p = c(p1,p2,p3)))
	  B1 <- B[1:1000,, drop = FALSE]
	  B2 <- B[1001:2023,, drop = FALSE]
	  B3 <- B[2024:3046,, drop = FALSE]

		k <- ncol(B)
		assignB <- attr(B, "assign")
		termlabelsB <- attr(terms(formula.p), "term.labels")
		coefnamesB <- colnames(B)
	}
	if(missing(s)){s <- matrix(1,q,k)}
	else{
		if(any(dim(s) != c(q,k))){stop("wrong size of 's'")}
		if(any(s != 0 & s != 1)){stop("'s' can only contain 0 and 1")}
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
		if(length(varsB) == 0){stop("the quantile function must depend on p")}
		if(length(constB) > 1){stop("remove multiple constant functions from 'formula.p'")}
		if(any(sB == 0 & mB == 0)){stop("remove zero functions from 'formula.p'")}
		sB[constB] <- B3[1,constB]; mB[constB] <- 0
	}
	else{
		sB <- rep.int(1, k); mB <- rep.int(0, k)
		if(intB){constB <- 1; varsB <- 2:k}
		else{constB <- integer(0); varsB <- 1:k}
	}

	if(all(s[, varsB] == 0)){stop("the quantile function must depend on p (wrong specification of 's')")}
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

	U <- list(X = X, y = y, z = z)
	X <- scale(X, center = mX, scale = sX)
	y <- (y - my)/(My - my)*10
	z <- z0 <- (z - my)/(My - my)*10
	z[z < min(y)] <- -Inf
	if(!use.slp){B3 <- scale(B3, center = mB, scale = sB)}

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
			}
		}
		MB <- colMeans(B3); MB[mB == 0] <- 0
		SB <- apply(B3,2,sd); SB[constB] <- 1
		B3 <- scale(B3, center = MB, scale = SB)
	}

	# Create a pre-evaluated basis (only used internally)

	p <- p3
	
	if(!use.slp){
	  bp <- B3
	  dp <- p[-1] - p[-1023]
	  b1p <- num.fun(dp,bp, "der")
	  Bp <- num.fun(dp,bp, "int")
	  BBp <- num.fun(dp,Bp, "int")
	  BB1 <- BBp[1023,]
	}
	else{
	  k <- attr(bfun, "k")
	  pp <- matrix(, 1023, k + 1)
	  pp[,1] <- 1; pp[,2] <- p
	  if(k > 1){for(j in 2:k){pp[,j + 1] <- pp[,j]*p}}
	  bp <- tcrossprod(pp, t(bfun$a))
	  b1p <- cbind(0, tcrossprod(pp[,1:k, drop = FALSE], t(bfun$a1[-1,-1, drop = FALSE])))
	  pp <- cbind(pp, pp[,k + 1]*p, pp[,k + 1]*p^2)
	  Bp <- tcrossprod(pp[,2:(k + 2)], t(bfun$A))
	  BBp <- tcrossprod(pp[,3:(k + 3)], t(bfun$AA))
	  BB1 <- colSums(bfun$AA)
	    
	  if(!intB){
	    bp <- bp[,-1, drop = FALSE]
	    b1p <- b1p[,-1, drop = FALSE]
	    Bp <- Bp[,-1, drop = FALSE]
	    BBp <- BBp[,-1, drop = FALSE]
	    BB1 <- BB1[-1]
	  }
	}
	BB1 <- matrix(rep(BB1, each = n), n)
	bpij <- NULL; for(i in 1:ncol(bp)){bpij <- cbind(bpij, bp*bp[,i])}

	internal.bfun <- list(p = p, bp = bp, b1p = b1p, Bp = Bp, BBp = BBp, BB1 = BB1, bpij = bpij)
        attr(internal.bfun, "pfun") <- approxfun(c(p[1], 0.5*(p[-1023] + p[-1])),p, method = "constant", rule = 2)

	# output. U = the original variables. V = the scaled/rotated variables.
	# stats.B, stats.X, stats.y = lists with the values use to scale/rotate

	stats.B <- list(m = mB, s = sB, M = MB, S = SB, rot = rotB, const = constB, vars = varsB,
		intercept = intB, term.labels = termlabelsB, assign = assignB, coef.names = coefnamesB)
	stats.X <- list(m = mX, s = sX, M = MX, S = SX, rot = rotX, const = constX, vars = varsX,
		intercept = intX, term.labels = termlabelsX, assign = assignX, coef.names = coefnamesX)
	stats.y <- list(m = my, M = My)

	V <- list(X = X, Xw = X*weights, y = y, z = z, d = d, weights = weights)
	if(type == "ctiqr"){V$z0 <- z0}
	list(mf = mf, U = U, V = V, stats.B = stats.B, stats.X = stats.X, stats.y = stats.y, 
		internal.bfun = internal.bfun, bfun = bfun, s = s, type = type)
}




#############################################################################################################
#############################################################################################################
#############################################################################################################


ctiqr.internal <- function(mf,cl, formula.p, tol = 1e-6, maxit, s, remove.qc){

	A <- check.in.iqr(mf, formula.p, s)
	V <- A$V; U <- A$U; s <- A$s; type <- A$type
	mf <- A$mf; n <- nrow(mf)
	S <- list(B = A$stats.B, X = A$stats.X, y = A$stats.y)
	attributes(A$bfun) <- c(attributes(A$bfun), S$B)
	bfun <- A$internal.bfun
	
	if(is.logical(remove.qc)){qc.rm <- remove.qc; qcc <- qc.control()}
	else{qc.rm <- TRUE; qcc <- remove.qc}

	if(missing(maxit)){maxit <- 10 + 10*sum(s)}
	else{maxit <- max(1, maxit)}
	
	
	q <- length(S$X$vars)
	if(type != "iqr" | q > 0){
	  Ty <- trans(V$z,V$y,V$d,V$weights, type)
	  yy <- Ty$f(V$y)
	  zz <- (if(type == "ctiqr") Ty$f(V$z) else V$z)
	}
	else{yy <- zz <- NULL}
	
	theta0 <- start.iqr(V$y, V$z, V$d, V$X, V$weights, bfun, 
	   df = max(5, min(15, round(n/30/(q + 1)))), yy, zz, s = s)
	
	Fit <- NULL
	fit.ok <- FALSE
	safeit <- 10
	try.count <- 0
	eeTol <- 1
	eps0 <- 0.05
	
	while(!fit.ok){
	  
	  try.count <- try.count + 1

		fit <- iqr.newton(theta0, V$y, V$z, V$d, V$X, V$Xw, 
			bfun, s, type, tol, maxit, safeit, eps0)
	
		if(max(abs(fit$ee)) > eeTol & try.count > 2)
		  {fit <- divide.et.impera(fit, V, bfun, s, type, tol, maxit, safeit, eps0)}

		if(fit.ok <- (fit$rank == ncol(fit$jacobian) & max(abs(fit$ee)) < eeTol)){
			covar <- try(cov.fun.iqr(fit$coefficients, V$y, V$z, V$d, V$X, V$Xw,
				V$weights, bfun, fit$p.star.y, fit$p.star.z, type, s = s), silent = TRUE)

			fit.ok <- (!inherits(covar, "try-error"))
		}
		
		covar.ok <- (if(fit.ok){!inherits(try(chol(covar$Q0), silent = TRUE), "try-error")} else FALSE)
		if(fit.ok & covar.ok){break}
		else if(fit.ok){Fit <- fit; Cov <- covar}
		
		if(try.count > 10 && !is.null(Fit)){break}
		if(try.count == 20){break}
		eeTol <- eeTol * 2 # I must be liberal. With discrete data, the objective function is nonsmooth.
		safeit <- safeit + 2
		eps0 <- eps0/2
	}

	if(!fit.ok && is.null(Fit)){stop("Unable to fit the model: 
		this can be due to severe misspecification, severe censoring and truncation,
		or overfitting (especially with discrete data)")}
	if(!covar.ok){warning("the estimated covariance matrix is deemed to be singular")}
	if(!fit.ok){fit <- Fit; covar <- Cov}
	if(!fit$converged){warning("the algorithm did not converge")}
	nit <- fit$n.it # I save it here, because if I remove quantile crossing, fit$n.it will change.

  # Quantile crossing

	if(qc.rm){
	  bfun$p.short <- bfun$p; bfun$bp.short <- bfun$bp; bfun$b1p.short <- bfun$b1p
	  pen0 <- qc.penalty(fit$coefficients, V$X, bfun, 1, H = FALSE)
	  if(pen0$pen == 0){warning("No detectable crossing, remove.qc is ignored"); qc.rm <- FALSE}
	}
  if(qc.rm){
	  if(type == "iqr"){loss0 <- iobjfun(fit$coefficients, V$y,V$X, V$weights, bfun, fit$p.star.y)}
	  else{loss0 <- iobjfun.ct(fit$coefficients, V$z,V$y,V$d,V$X,V$weights, bfun, fit$py, fit$pz, type)}
    
    if(lambda.in <- (!is.null(qcc$lambda))){lambda <- qcc$lambda; qcc$maxTry <- 1}
	  else{lambda <- 1/abs(pen0$pen)*0.25*loss0} # Initial total penalization = 0.25*loss
	  
	  # I add information to "fit". The goal is: if fixqc fails, at least I can select
	  # the initial fit as my "best" fit. Also, to use fixqc, I need b''(p).
	  fit$covar <- covar
	  fit$covar.ok <- covar.ok
	  fit$pen <- pen0 # before qr.cm, fit$pen was 0 (as lambda was 0)
	  dp <- bfun$p[-1] - bfun$p[-1023]
	  bfun$b2p <- num.fun(dp,bfun$b1p, "der")

	  if(!lambda.in){
	    r <- c(50, 100, 250, 500, 1023)
	    bycross <- c(2,1,1,1,1)
	    tol <- min(tol, 1e-6) # "tol" should not affect the ability to remove crossing.
	  }
	  else{r = 1023; bycross <- 1; tol <- 1e-10}
	  
	  
	  count <- 0
	  lambda0 <- lambda
	  for(i in 1:length(r)){

	    eps0 <- (if(i == 1) 0.1 else fit$eps*32)
	    
	    pcross <- which(pen0$pcross) # quantiles at which crossing occurs
	    fit <- fixqc(fit, V, bfun, s, type, tol, maxit = maxit, safeit = safeit, eps0 = eps0,
	      lambda = lambda, r = r[i], maxTry = qcc$maxTry, trace = qcc$trace, 
	      count = count, pcross = pcross[seq.int(1,length(pcross), bycross[i])])
	    # remark: I restart from the last used eps0. Large eps ---> very wrong parameters ---> 
	    # a lot of obs with crossing ---> qc.penalty becomes much much slower!
	    
	    if(!fit$warn & fit$done){ # If you believe you did it, let's check if you really did it.
	      pen0 <- qc.penalty(fit$coefficients, V$X, bfun, 1, H = FALSE)
	      if(pen0$pen == 0){if(qcc$trace){cat("###################", "\n"); cat("Success.", "\n")}; break}
	      else{fit$done <- FALSE}
	    }

	    count <- fit$count
	    if(count == qcc$maxTry){
	        if(!lambda.in){warning("Quantile crossing could not be removed entirely, please increase 'maxTry'")}
	        else{warning("Quantile crossing could not be removed entirely at the provided value of lambda")}
	        break
	    }
	    
	    if(i == length(r) && !fit$done){
	      if(!lambda.in){warning("Quantile crossing could not be removed entirely")}
	      else{warning("Quantile crossing could not be removed entirely at the provided value of lambda")}
	    }

	    lambda0 <- fit$lambda # I restart from the last lambda, increasing it a bit (see below)
	    lambda <- lambda*1.5
	    tol <- tol/2 # failure may be caused by too large tol
	  }
	  if(!fit$covar.ok){warning("After remove.qc: the estimated covariance matrix is deemed to be singular")}
	  if(!fit$converged){warning("After remove.qc: the algorithm did not converge")}
	   
    covar <- fit$covar
  }

	# minimized loss function

	if(type == "iqr"){obj <- iobjfun(fit$coef, V$y,V$X,V$weights, bfun, fit$p.star.y)}
	else{obj <- iobjfun.ct(fit$coef, V$z,V$y,V$d,V$X,V$weights, bfun, fit$py, fit$pz, type)}

	obj <- obj*(S$y$M - S$y$m)/10
	attr(obj, "df") <- sum(s)
	if(type != "iqr"){attr(obj, "message") <- "Not the function being minimized"}
	fit$obj.function <- obj

	# fitted CDFs (for internal use, precision ~ 0.001)

	CDFs <- data.frame(CDF.y = fit$py, 
	CDF.z = (if(type == "ctiqr") fit$pz else NA))
	attr(CDFs, "km") <- km(CDFs$CDF.z, CDFs$CDF.y, V$d, V$weights, type)

	# output
  
	attr(mf, "assign") <- S$X$assign
	attr(mf, "stats") <- S
	attr(mf, "CDFs") <- CDFs
	attr(mf, "all.vars") <- V
	attr(mf, "all.vars.unscaled") <- U
	attr(mf, "Q0") <- covar$Q0
	attr(mf, "internal.bfun") <- bfun
	attr(mf, "bfun") <- A$bfun
	attr(mf, "theta") <- fit$coefficients
	attr(mf, "type") <- type

	out <- check.out(fit$coefficients, S, covar = covar$Q)
	fit <- list(coefficients = out$theta, call = cl, 
		converged = fit$converged, n.it = nit,
		obj.function = fit$obj.function,
		covar = out$covar, mf = mf, s = s)

	jnames <- c(sapply(attr(A$bfun, "coef.names"), 
		function(x,y){paste(x,y, sep = ":")}, y = S$X$coef.names))
	dimnames(fit$covar) <- list(jnames, jnames)
	dimnames(fit$coefficients) <- dimnames(fit$s) <- list(S$X$coef.names, S$B$coef.names)

	# CDF and PDF, precision ~ 1e-6. I avoid CDF < 2e-6 and CDF > 1 - 2e-6, because
	# these two are the extreme of the grid used internally. What happens beyond these
	# two values is not under control, and for example there could be crossing.

	fit$CDF <- pmin(1 - 2e-6, pmax(2e-6, p.bisec(fit$coef,U$y,U$X,A$bfun)))
	b1 <- apply_bfun(A$bfun, fit$CDF, "b1fun")
	fit$PDF <- 1/c(rowSums((U$X%*%fit$coef)*b1))
	if(any(fit$PDF < 0)){warning("Quantile crossing detected (PDF < 0)")}
	# fit$PDF[attr(fit$CDF, "out")] <- 0 # removed in v3.0
	attributes(fit$CDF) <- attributes(fit$PDF) <- list(names = rownames(mf))
	

	# finish

	class(fit) <- "iqr"
	fit
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


iobjfun <- function(theta, y,X,weights, bfun, p.star){
	s <- NULL
	B <- bfun$Bp[p.star,, drop = FALSE]
	py <- bfun$p[p.star]
	for(j in 1:ncol(B)){s <- cbind(s, bfun$BB1[1,j] - B[,j])}
	sum(y*weights*(py - 0.5)) + sum(((X*weights)%*%theta)*s)
}

iobjfun.ct <- function(theta, z,y,d,X,weights, bfun, py, pz, type){
 
  trunc <- (type == "ctiqr")
  n <- nrow(X)
  
  ###########################################################
  
  # I use a shorter grid
  r <- 250
  psel <- round(seq.int(1,1023, length.out = r))
  p <- bfun$p[psel]
  bp <- bfun$bp[psel,, drop = FALSE]
  dp <- p[-1] - p[-r]; dp <- c(dp[1], dp)
  
  ###########################################################
  
  beta <- tcrossprod(theta, bp)
  eta <- tcrossprod(X, t(beta))
  
  deltay <- y - eta
  omegay <- (deltay <= 0)
  Sy <- 1 - matrix(pmin(py, 1-1e-6), n, r)
  
  deltaz <- (if(trunc) z - eta else -Inf)
  omegaz <- (if(trunc) deltaz <= 0 else 1)
  Sz <- 1 - (if(trunc) matrix(pmin(pz, 1-1e-6), n, r) else 0)
  
  ###########################################################
  
  p <- t(matrix(t(p),r,n))
  dp <- t(matrix(t(dp),r,n))
  pbar <- 1 - p
  
  omega.hat <- omegay*(1 - (1 - d)*pbar/Sy) - omegaz*(1 - pbar/Sz) # this is actually (omega.hat - p)
  loss <- -weights*(omega.hat*deltay)
  loss <- .rowSums(loss*dp, n,r)
  sum(loss)
}

#############################################################################################################
#############################################################################################################
#############################################################################################################

iqr.ee <- function(theta, y,z,d,X,Xw, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE, lambda = 0){

	k <- ncol(theta)
	n <- length(y)
	BB1 <- bfun$BB1
	if(missing(G)){
		B <- bfun$Bp[p.star.y,, drop = FALSE]
		S1 <- BB1 - B
		pen <- qc.penalty(theta, X, bfun, lambda, H = FALSE)

		if(!i){g <- c(crossprod(Xw,S1)) - lambda*pen$gradient}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B <- G$B; g <- G$g; pen <- G$pen}

	if(J){
		b1 <- bfun$b1p[p.star.y,, drop = FALSE]
		bij <- bfun$bpij[p.star.y,, drop = FALSE]

		A1 <- 1/c(.rowSums(tcrossprod(X, t(theta))*b1, n,k))
		A1 <- pmax0(A1)
		A1[attr(p.star.y, "out")] <- 0
		Xw <- Xw*A1

		J <- NULL
		count <- 0
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){
				count <- count + 1
				h.temp <- cbind(h.temp, crossprod(Xw, X*bij[,count]))
			}
			J <- rbind(J, h.temp)
		}
		
    pen <- qc.penalty(theta, X, bfun, lambda, pen = pen, H = TRUE)
		J <- J - lambda*pen$hessian
	}
	
	list(g = g, J = J, B = B, pen = pen)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


ciqr.ee <- function(theta, y,z,d,X,Xw, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE, lambda = 0){

	k <- ncol(theta)
	n <- length(y)
	BB1 <- bfun$BB1

	if(missing(G)){
		B <- bfun$Bp[p.star.y,, drop = FALSE]
		BB <- bfun$BBp[p.star.y,, drop = FALSE]
		py <- bfun$p[p.star.y]

		a <- (1 - d)/(1 - py); a[attr(p.star.y, "out.r")] <- 0
		S1 <- a*(BB - BB1)
		a <- d; a[attr(p.star.y, "out.r")] <- 1
		S1 <- BB1 - a*B + S1
		pen <- qc.penalty(theta, X, bfun, lambda, H = FALSE)

		if(!i){g <- c(crossprod(Xw,S1)) - lambda*pen$gradient}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B <- G$B; BB <- G$BB; g <- G$g; py <- G$py; pen <- G$pen}

	if(J){
		b <- bfun$bp[p.star.y,, drop = FALSE]
		b1 <- bfun$b1p[p.star.y,, drop = FALSE]

		A1 <- 1/c(.rowSums(tcrossprod(X, t(theta))*b1, n,k))
		A1 <- pmax0(A1)
		A1[attr(p.star.y, "out")] <- 0
		a <- 1 - py
		a[attr(p.star.y, "out.r")] <- 1 
		# the value 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1)
		A2 <- d*b + (1 - d)/(a^2)*(BB1 - B*a - BB)
		Xw <- Xw*A1

		J <- NULL
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){h.temp <- cbind(h.temp, crossprod(Xw, X*(b[,i2]*A2[,i1])))}
			J <- rbind(J, h.temp)
		}
		pen <- qc.penalty(theta, X, bfun, lambda, pen = pen, H = TRUE)
		J <- J - lambda*pen$hessian
	}
	list(g = g, J = J, B = B, BB = BB, py = py, pen = pen)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


ctiqr.ee <- function(theta, y,z,d,X,Xw, bfun, 
	p.star.y, p.star.z, J = TRUE, G, i = FALSE, lambda = 0){

	k <- ncol(theta)
	n <- length(y)
	BB1 <- bfun$BB1
	if(missing(G)){
		B.y <- bfun$Bp[p.star.y,, drop = FALSE]	
		BB.y <- bfun$BBp[p.star.y,, drop = FALSE]
		B.z <- bfun$Bp[p.star.z,, drop = FALSE]
		BB.z <- bfun$BBp[p.star.z,, drop = FALSE]
		py <- bfun$p[p.star.y]
		pz <- bfun$p[p.star.z]
	
		out.y <- attr(p.star.y, "out.r")
		out.z <- attr(p.star.z, "out.r")
		a <- d

		S1.y <- (1 - d)/(1 - py)*(BB.y - BB1)
		S1.z <- 1/(1 - pz)*(BB.z - BB1)
		S1.y[out.y,] <- 0; a[out.y] <- 1
		S1.y[out.z,] <- S1.z[out.z,] <- a[out.z] <- 0

		S1 <- -a*B.y + S1.y - S1.z
		pen <- qc.penalty(theta, X, bfun, lambda, H = FALSE)
		
		if(!i){g <- c(crossprod(Xw,S1)) - lambda*pen$gradient}
		else{g <- NULL; for(h in 1:k){g <- cbind(g,X*S1[,h])}}
	}
	else{B.y <- G$B.y; BB.y <- G$BB.y; B.z <- G$B.z; BB.z <- G$BB.z; g <- G$g; py <- G$py; pz <- G$pz; pen <- G$pen}

	if(J){
		b.y <- bfun$bp[p.star.y,, drop = FALSE]
		b1.y <- bfun$b1p[p.star.y,, drop = FALSE]
		b.z <- bfun$bp[p.star.z,, drop = FALSE]
		b1.z <- bfun$b1p[p.star.z,, drop = FALSE]

		Xtheta <- tcrossprod(X, t(theta))

		A1.y <- 1/c(.rowSums(Xtheta*b1.y, n,k))
		A1.z <- 1/c(.rowSums(Xtheta*b1.z, n,k))
		A1.y <- pmax0(A1.y)
		A1.z <- pmax0(A1.z)
		A1.y[attr(p.star.y, "out")] <- 0
		A1.z[attr(p.star.z, "out")] <- 0

		ay <- 1 - py; az <- 1 - pz
		ay[attr(p.star.y, "out.r")] <- az[attr(p.star.z, "out.r")] <- 1 
		# the value 1 is arbitrary, only to avoid NAs (will be canceled by the zeroes in A1.y and A1.z)

		A2.y <- d*b.y - (1 - d)/(ay^2)*(B.y*ay + BB.y - BB1)
		A2.z <- 1/(az^2)*(B.z*az + BB.z - BB1)
		H.y <- A2.y*A1.y 
		H.z <- A2.z*A1.z

		J <- NULL
		for(i1 in 1:k){
			h.temp <- NULL
			for(i2 in 1:k){h.temp <- cbind(h.temp, crossprod(Xw, X*(b.y[,i2]*H.y[,i1] + b.z[,i2]*H.z[,i1])))}
			J <- rbind(J, h.temp)
		}
		pen <- qc.penalty(theta, X, bfun, lambda, pen = pen, H = TRUE)
		J <- J - lambda*pen$hessian
	}
	list(g = g, J = J, B.y = B.y, BB.y = BB.y, B.z = B.z, BB.z = BB.z, py = py, pz = pz, pen = pen)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


start.iqr <- function(y,z,d, x, weights, bfun, df, yy, zz, s){

	if(is.null(yy)){p.star <- (rank(y, ties.method = "first") - 0.5)/length(y)}
	else{
	  pch.fit.ct <- getFromNamespace("pch.fit.ct", ns = "pch")
	  predF.pch <- getFromNamespace("predF.pch", ns = "pch")
	  m0 <- suppressWarnings(pch.fit.ct(z = zz, y = yy, d = d, 
		x = cbind(1,x), w = weights, breaks = df))
	  p.star <- 1 - predF.pch(m0)[,3]
	}

	pfun <- attr(bfun, "pfun")
	p.star <- pfun(p.star)  
	b.star <- bfun$bp[match(p.star, bfun$p),]
	X <- model.matrix(~ -1 + x:b.star)
	X <- X[, c(s) == 1, drop = FALSE]
	start.ok <- FALSE; restol <- 4
	while(!start.ok){
	  
		m <- lm.wfit(X, y, weights)
		res <- m$residuals
		start.ok <- all(w <- (abs(res)/sd(res) < restol))

		if(start.ok | sum(w) < 0.5*length(y)){break}
		X <- X[w,, drop = FALSE]
		y <- y[w]
		weights <- weights[w]
		restol <- restol + 1
	}
	out <- rep.int(0, length(s))
	out[c(s) == 1] <- m$coef
	out <- matrix(out, ncol(x))
	out[is.na(out)] <- 0
	out
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


# Note: if s has some zeroes, the covariance matrix Q will contain some zero-columns and rows,
# while the gradient and jacobian will just omit the parameters that are not estimated

cov.fun.iqr <- function(theta, y,z,d,X,Xw, weights, bfun, 
		p.star.y, p.star.z, type, s){

	if(type == "iqr"){ee <- iqr.ee}
	else if(type == "ciqr"){ee <- ciqr.ee}
	else{ee <- ctiqr.ee}
	s <- c(s == 1)
	
	G.i <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = TRUE, i = TRUE)
	s.i <- G.i$g[,s, drop = FALSE]
	G <- G.i$J[s,s, drop = FALSE]
	
	# to avoid singularities/nondifferentiability (added in v 3.0)
	s.i <- t(t(s.i) - colMeans(s.i))
	G <- G + diag(0.01, ncol(G))

	Omega <- chol2inv(chol(t(s.i*weights)%*%s.i + diag(0.01, ncol(s.i)))) # extra diag added in v 3.0
	Q0 <- t(G)%*%Omega%*%G + diag(0.01, ncol(G)) # extra diag added in v 3.0
	Q <- chol2inv(chol(Q0))

	# I reduce correlation between estimates.
	# If it is too high, something may go wrong: I bound it at 0.999.

	if((cc <- ncol(Q)) > 1){
	  for(j1 in 1:(cc - 1)){
	    for(j2 in (j1 + 1):cc){
	      s12 <- sign(Q[j1,j2])
	      Q[j1,j2] <- Q[j2,j1] <- s12*pmin(abs(Q[j1,j2]), 0.999*sqrt(Q[j1,j1]*Q[j2,j2])) 
	    }
	  }
	}

	# Finish
	
	U <- matrix(0, length(s), length(s))
	U[s,s] <- Q
	list(Q = U, Q0 = Q0, jacobian = G, ee = colSums(s.i*weights), Omega = Omega, s.i = s.i)
}






#############################################################################################################
#############################################################################################################
#############################################################################################################

iqr.newton <- function(theta, y,z,d,X,Xw, bfun, s, type, tol, maxit, safeit, eps0, lambda = 0){
 
  if(type == "iqr"){ee <- iqr.ee}
  else if(type == "ciqr"){ee <- ciqr.ee}
  else{ee <- ctiqr.ee}
  
  q <- nrow(theta)
  k <- ncol(theta)
  s <- c(s == 1)
  
  p.star.y <- p.bisec.internal(theta, y,X, bfun$bp)
  if(type == "ctiqr"){p.star.z <- p.bisec.internal(theta, z,X, bfun$bp)}
  G <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, lambda = lambda)
  pen <- G$pen
  
  g <- G$g[s]
  conv <- FALSE
  eps <- eps0

  # Preliminary safe iterations, only g is used

  for(i in 1:safeit){

    if(conv | max(abs(g)) < tol){break}
    u <- rep.int(0, q*k)
    u[s] <- g
    delta <- matrix(u, q,k)	
    delta[is.na(delta)] <- 0
    cond <- FALSE

    while(!cond){
      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p.star.y <- p.bisec.internal(new.theta, y,X, bfun$bp)
      if(type == "ctiqr"){p.star.z <- p.bisec.internal(new.theta, z,X, bfun$bp)}
      G1 <- ee(new.theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, lambda = lambda)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }

    if(conv){break}
    g <- g1
    G <- G1
    theta <- new.theta
    eps <- min(eps*2,0.1)
  }

  # Newton-Raphson

  alg <- "nr"
  conv <- FALSE
  eps <- eps*10
  h <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = TRUE, G = G, lambda = lambda)$J[s,s, drop = FALSE]
  h <- h + diag(0.0001, nrow(h)) # added in version 3.0
 
  for(i in 1:maxit){

    if(conv | max(abs(g)) < tol){break}
    
    ####
    
    if(type == "iqr"){
      H1 <- try(chol(h), silent = TRUE)
      err <- inherits(H1, "try-error")
    }
    else{
      H1 <- qr(h)
      r <- H1$rank
      err <- (r != ncol(h))
    }
    if(!err){
      if(alg == "gs"){alg <- "nr"; eps <- 1}
      delta <- (if(type == "iqr") chol2inv(H1)%*%g else qr.solve(H1)%*%g)
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
      new.theta <- theta - delta*eps
      if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
      p.star.y <- p.bisec.internal(new.theta, y,X, bfun$bp)
      if(type == "ctiqr"){p.star.z <- p.bisec.internal(new.theta, z,X, bfun$bp)}
      G1 <- ee(new.theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = FALSE, lambda = lambda)
      g1 <- G1$g[s]
      cond <- (sum(g1^2) < sum(g^2))
      eps <- eps*0.5
    }
    
    if(conv){break}

    g <- g1
    G <- G1
    theta <- new.theta
    h <- ee(theta, y,z,d,X,Xw, bfun, p.star.y, p.star.z, J = TRUE, G = G, lambda = lambda)$J[s,s, drop = FALSE]
    if(i > 1){eps <- min(eps*10,1)}
    else{eps <- min(eps*10,0.1)}
  }
  
  p.star.y <- p.bisec.internal(theta, y,X, bfun$bp)
  py <- bfun$p[p.star.y]
  if(type == "ctiqr"){
    p.star.z <- p.bisec.internal(theta, z,X, bfun$bp)
    pz <- bfun$p[p.star.z]
    pz <- pmin(pz, py - 1e-8)
    pz[p.star.z == 1] <- 0
  }
  else{p.star.z <- pz <- NULL}

  list(coefficients = matrix(theta, q, k),
       converged = (i < maxit), n.it = i, eps = eps,
       p.star.y = p.star.y, p.star.z = p.star.z, py = py, pz = pz,
       ee = g, jacobian = h, pen = G$pen,
       rank = (alg == "nr")*sum(s))
}


# Fix the max ee. If fail, try the second largest, and so on.
divide.et.impera <- function(fit, V, bfun, s, type, tol, maxit, safeit, eps0, lambda = 0){
  
  if(sum(s) == 1){fit$failes <- FALSE; return(fit)}
  theta0 <- theta <- fit$coef
  ee <- s; ee[s == 1] <- abs(fit$ee)
  
  failed <- TRUE
  while(failed){

    i <- maxind(ee)

    sbis <- s*(ee > 0); sbis[-i[1],] <- 0
    fit <- iqr.newton(theta, V$y, V$z, V$d, V$X, V$Xw, 
      bfun, sbis, type, tol = tol, maxit = 2*sum(sbis), safeit = 5, eps0 = eps0, lambda = lambda)
    theta <- fit$coef

    sbis <- s*(ee > 0); sbis[,-i[2]] <- 0
    fit <- iqr.newton(theta, V$y, V$z, V$d, V$X, V$Xw, 
      bfun, sbis, type, tol = tol, maxit = 2*sum(sbis), safeit = 5, eps0 = eps0, lambda = lambda)
    theta <- fit$coef
  
    failed <- (all(theta == theta0))
    ee[i[1],i[2]] <- 0 # If I failed, I just give up this coefficient
    if(all(ee == 0)){break}
  }

  
  if(failed){fit$failed <- TRUE; return(fit)}
  
  fit <- iqr.newton(theta, V$y, V$z, V$d, V$X, V$Xw, 
    bfun, s, type, tol, maxit, safeit, eps0, lambda)
  fit$failed <- FALSE
  fit
}


