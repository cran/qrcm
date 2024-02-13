#' @export
diagnose.qc <- function(obj){
  
  if(!inherits(obj, "iqr")){stop("this function can only be applied to 'iqr' objects")}
  fittype <- attr(attr(obj$mf, "type"), "fittype")
  if(fittype == "iciqr"){ # I take pmin(PDF$left, PDF$right) and select CDF accordingly
    w <- (obj$PDF[,1] > obj$PDF[,2]) + 1
    w <- cbind(1:length(w), w)
    obj$PDF <- obj$PDF[w]
    obj$CDF <- obj$CDF[w]
  }
  
  bfun <- attr(obj$mf, "internal.bfun")
  X <- attr(obj$mf, "all.vars")$X
  theta <- attr(obj$mf, "theta")
  n <- nrow(X); r <- length(bfun$p)

  ########################################################################################
  
  Q1 <- tcrossprod(tcrossprod(X, t(theta)), bfun$b1p)
  sel <- which(.rowSums((sQ1 <- (Q1 < 0)),n,r) > 0)
  nsel <- length(sel)
  
  ########################################################################################

  if(nsel == 0){
    out <- list(
      qc = data.frame(qc.local = rep.int(FALSE,n), qc.global = rep.int(FALSE,n)),
      qc.local = 0, qc.global = 0, pcross = NULL, crossIndex = 0
    )
    class(out) <- "qc.iqr"
    return(out)
  }
  
  sQ1 <- sQ1[sel,, drop = FALSE]
  pleft <- pright <- sQ1
  for(j in 2:r){pleft[,j] <- pleft[,j]*(1 - sQ1[,j-1])}
  for(j in (r - 1):1){pright[,j] <- pright[,j]*(1 - sQ1[,j+1])}
  single_p <- (pleft & pright) # where Q'(p) < 0 at a single value of p.
  
  ########################################################################################

  # Local and global quantile crossing
  
  qc.global <- rep.int(FALSE, n); qc.global[sel] <- apply(sQ1,1,any)
  qc.local <- (obj$PDF < 0)
  
  undetected.crossing <- (qc.local & !qc.global)
  qc.global[undetected.crossing] <- TRUE
  
  ########################################################################################
  
  # A couple of measures of crossing.
  
  # (1) Quantiles at which there is quantile crossing

  if(any(sQ1) | any(undetected.crossing)){
    pcross <- rep.int(1:r, .colSums(sQ1, nsel,r))
    pcross <- c(bfun$p[pcross], obj$CDF[undetected.crossing])
    pcross <- cut(pcross, c(0,0.001,0.01,0.05,0.10,0.25,0.5,0.75,0.90,0.95,0.99,0.999,1), include.lowest = TRUE)
    pcross <- table(pcross)/length(pcross)
    pcross <- cbind("%" = pcross*100)
  }
  else{pcross <- NULL}
  
  # (2) A scale-invariant indicator of how much quantile crossing there is:
  # average proportion of the (0,1) line where quantiles cross

  deltap <- pright - pleft # (p+) - (p-)
  crossIndex <- tcrossprod(deltap, t(bfun$p))
  crossIndex <- (sum(crossIndex) + 1e-6*sum(single_p) + 1e-6*sum(undetected.crossing))/n 
  # For each single-p crossing or undetected crossing, I assume the interval has length 1e-6 
  # [which is less than the shortest distance between any two values of bfun$p]
  
  ########################################################################################
  # output ###############################################################################
  ########################################################################################

  qc <- data.frame(qc.local = qc.local, qc.global = qc.global)
  rownames(qc) <- rownames(obj$mf)

  out <- list(qc = qc, qc.local = sum(qc$qc.local), qc.global = sum(qc$qc.global), 
    pcross = pcross, crossIndex = crossIndex)
  class(out) <- "qc.iqr"
  out
}


#' @export
print.qc.iqr <- function(x, ...){
  
  n <- nrow(x$qc)
  
  n.qc <- sum(x$qc.local)
  p.qc <- round(n.qc/n*100,3)
  cat("Local quantile crossing [f(y | x, theta) < 0]:", "\n")
  cat("n. of observations: ", n.qc, " (", p.qc, "%)", "\n", sep = "")
  
  cat("\n")
  
  n.qc <- sum(x$qc.global)
  p.qc <- round(n.qc/n*100,3)
  cat("Global quantile crossing [Q'(p | x, theta) < 0 at some p]:", "\n")
  cat("n. of observations: ", n.qc, " (", p.qc, "%)", "\n", sep = "")
  
  cat("\n")
  
  if(!is.null(x$pcross)){
    
    cat("Distribution of p at which global quantile crossing occurs:", "\n")
    print(round(x$pcross,2))
    cat("\n")
    
    cat("Cross Index", "\n")
    cat(signif(x$crossIndex,4), "\n\n")
  }
  else{cat("Cross Index", "\n"); cat(0, "\n\n")}
}


qc.penalty <- function(theta, X, bfun, lambda, pen, H){
  
  if(lambda == 0){return(list(pen = 0, gradient = 0, hessian = 0))}

  q <- nrow(theta); k <- ncol(theta); N <- n <- nrow(X)
  p <- bfun$p.short; r <- length(p)
  bp <- bfun$bp.short; b1p <- bfun$b1p.short; b2p <- bfun$b2p.short
  
  if(missing(pen)){
    Xtheta <- tcrossprod(X, t(theta))
    Q1 <- tcrossprod(Xtheta, b1p)
    sQ1 <- (Q1 < 0)

    sel <- which(.rowSums(sQ1,n,r) > 0)
    sQ1 <- sQ1[sel,, drop = FALSE]
    X <- X[sel,, drop = FALSE]
    Xtheta <- Xtheta[sel,, drop = FALSE]
    n <- length(sel)
  
    pleft <- pright <- sQ1
    for(j in 2:r){pleft[,j] <- pleft[,j]*(1 - sQ1[,j-1])}
    for(j in (r - 1):1){pright[,j] <- pright[,j]*(1 - sQ1[,j+1])}

    deltap <- pright - pleft # (p+) - (p-)
    deltab <- tcrossprod(deltap, t(bp)) # b(p+) - b(p-)
    
    if(any(single_p <- (pleft & pright))){
      deltab <- deltab + 1e-6*tcrossprod(single_p, t(b1p))
    }
    
    crossIndex <- tcrossprod(deltap, t(p)) # proportion of the (0,1) line where quantile cross
    crossIndex <- (sum(crossIndex) + 1e-6*sum(single_p))/N # a scale-invariant indicator of how much quantile crossing there is.
    pen <- sum(.rowSums(Xtheta*deltab, n,k)) # penalization
    g <- c(crossprod(X, deltab)) # first derivative of the penalization
  }
  else{
    pleft <- pen$pleft; pright <- pen$pright; sQ1 <- pen$sQ1
    crossIndex <- pen$crossIndex; pcross <- pen$pcross
    n <- pen$n; X <- pen$X; Xtheta <- pen$Xtheta
    g <- pen$g; pen <- pen$pen
  }
  
  #########################################################################
  
  if(H){
    pleft[,1] <- pright[,r] <- 0 # points where Q' is not actually zero do not enter the hessian
    b1left <- tcrossprod(pleft, t(b1p)) # b'(pleft)
    b1right <- tcrossprod(pright, t(b1p)) # b'(pright)
    b2left <- tcrossprod(pleft, t(b2p)) # b''(pleft)
    b2right <- tcrossprod(pright, t(b2p)) # b''(pright)
    Q2left <- .rowSums(Xtheta*b2left, n,k) # Q''(pleft)
    Q2right <- .rowSums(Xtheta*b2right, n,k) # Q''(pright)

    H <- NULL
    count <- 0

    Xleft <- X/Q2left; Xleft[is.na(Xleft) | !is.finite(Xleft)] <- 0
    Xright <- X/Q2right; Xright[is.na(Xright) | !is.finite(Xright)] <- 0

    for(i1 in 1:k){
      h.left <- h.right <- NULL
      for(i2 in 1:k){
        h.left <- cbind(h.left, crossprod(Xleft, X*b1left[,i1]*b1left[,i2]))
        h.right <- cbind(h.right, crossprod(Xright, X*b1right[,i1]*b1right[,i2]))
      }
      H <- rbind(H, h.left - h.right)
    }
  }
  
  #########################################################################
  
  list(pen = pen, gradient = g, hessian = H, 
       pleft = pleft, pright = pright, sQ1 = sQ1,
       crossIndex = crossIndex, pcross = (.colSums(sQ1,n,r) > 0),
       n = n, X = X, Xtheta = Xtheta)
}

# The argument "pcross" includes the indexes of the values of p
# at which quantile crossing was detected by a previous call to qc.penalty.
# The goal is to use a grid which is more "dense" around the critical values of p.
fixqc <- function(fit, V, bfun, s, type, tol, maxit, safeit, eps0, 
  lambda, r, maxTry, trace, 
  count, pcross = NULL){
  
  # I create a shorter grid to estimate the penalty term

  psel <- c(round(seq.int(1,1023, length.out = r)), pcross)
  psel <- unique(sort(psel))
  bfun$p.short <- bfun$p[psel]
  bfun$bp.short <- bfun$bp[psel,, drop = FALSE]
  bfun$b1p.short <- bfun$b1p[psel,, drop = FALSE]
  bfun$b2p.short <- bfun$b2p[psel,, drop = FALSE]
  
  # Initial value of the penalization. If zero, I need a longer grid.
  theta0 <- fit$coefficients
  pen0 <- qc.penalty(theta0, V$X, bfun, 1, H = FALSE)
  if(pen0$pen == 0){
    fit$warn <- TRUE
    fit$done <- FALSE
    return(fit)
  }
  
  ##################################################################

  fit.ok <- qc.ok <- covar.ok <- FALSE
  lambdaMin <- 0; lambdaMax <- 2*lambda
  lambda_not_too_large <- TRUE; maxLambda_with_fit_ok <- 0
  
  lambda_fail <- theta_fail <- 0 
  # lambda_fail counts how many times in total I had a lambda that caused fit.ok = FALSE.
  # theta_fail counts how many times IN A ROW I could not move theta.
  
  Fit <- list(fit); Covar <- list(fit$covar); Covar.ok <- fit$covar.ok
  CrossIndex <- fit$pen$crossIndex; Lambda <- 1e-4 
  # The first lambda was actually zero, but I say 1e-4: if I fail to move,
  # I don't want to return a lambda = 0, because then the fit does not
  # create the "pen" object.
  
  
  done <- FALSE
  while(!(fit.ok & covar.ok & qc.ok)){

    count <- count + 1

    fit <- iqr.newton(theta0, V$y, V$z, V$d, V$X, V$Xw, 
        bfun, s, type, tol, maxit, safeit, eps0, lambda = lambda)
    theta_fail <- (if(all(fit$coefficients == theta0)) theta_fail + 1 else 0)

    # Check n.1: theta did not move, rank-deficient fit.
    if(theta_fail > 2 | (fit$rank != ncol(fit$jacobian) & lambda_fail > 5)){
        fit <- divide.et.impera(fit, V, bfun, s, type, tol, maxit, safeit, eps0, lambda)
        theta_fail <- 0
    }

    # Check n.2: estimating equation. 
    eeTol <- 1
    while(max(abs(fit$ee)) > eeTol){
      fit <- divide.et.impera(fit, V, bfun, s, type, tol, maxit, safeit, eps0, lambda)
      if(fit$failed){break} # actually, failing completely is very rare
      eeTol <- eeTol*1.5
    }

    # IMPORTANT: at convergence, the ee is **NOT** expected to be zero, being an unpenalized
    # estimating equation (if pen = 0, its gradient is also zero). This means that
    # sometimes I may have a "bad" ee which is actually ok. The above check (ee < 1),
    # however, will save the fit a lot of times. Sometimes, it will do nothing and just waste time.

    ##################################################################

    qc.ok <- (fit$pen$crossIndex == 0)
    fit.ok <- (fit$rank == ncol(fit$jacobian))

    if(fit.ok){
      covar <- try(cov.fun.iqr(fit$coefficients, V$y, V$z, V$d, V$X, V$Xw,
        V$weights, bfun, fit$p.star.y, fit$p.star.z, type, s = s), silent = FALSE)
      fit.ok <- !inherits(covar, "try-error")
    }
    
    covar.ok <- (if(fit.ok){!inherits(try(chol(covar$Q0), silent = TRUE), "try-error")} else FALSE)

    # Did I finish?

    if(trace){
      est.cI <- signif(fit$pen$crossIndex,6)
      if(est.cI == 0){est.cI <- signif(1e-6/nrow(V$X),6); est.cI <- paste("<", est.cI)}
      # if crossIndex is zero here, it does not mean that it is zero if I evaluate it more accurately.
      cat("###################", "\n")
      cat("Iteration", count, "of", maxTry, "\n")
      cat("lambda:", signif(lambda,4), "\n")
      cat("CrossIndex:", est.cI, "\n")
      cat("Full rank fit:", fit.ok & covar.ok, "\n")
      cat("\n")
    }

    if(fit.ok & covar.ok & qc.ok){
      done <- TRUE
      break
    }
    if(fit.ok){
      Fit[[length(Fit) + 1]] <- fit
      Covar[[length(Covar) + 1]] <- covar
      Covar.ok <- c(Covar.ok, covar.ok)
      CrossIndex <- c(CrossIndex, fit$pen$crossIndex)
      Lambda <- c(Lambda, lambda)
    }
    if(count == maxTry){break}
    if(!is.null(fit$failed) && fit$failed){break}
    
    # Assuming I did not finish, what do I do?
    
    if(fit.ok & covar.ok){
      theta0 <- fit$coefficients
      maxLambda_with_fit_ok <- max(maxLambda_with_fit_ok, lambda)
    }
    else{lambda_not_too_large <- FALSE}
    lambda_fail <- lambda_fail + !fit.ok

    if(!qc.ok){
      # if !fit.ok, for 5 times I still allow lambda to increase. Sometimes fit.ok becomes TRUE at larger lambdas.
      if(lambda_not_too_large | lambda_fail <= 5){lambdaMin <- lambda; lambdaMax <- 2*lambdaMax}
      else{
        lambdaMin <- maxLambda_with_fit_ok
        lambdaMax <- (if(fit.ok & covar.ok) 2*maxLambda_with_fit_ok else lambda)
      }
    }
    else{lambdaMax <- lambda}
    lambda <- (lambdaMin + lambdaMax)/2
  }
  
  if(!done){
    if(any(Covar.ok)){
      sel <- which(Covar.ok)
      if(length(sel) > 1){sel <- sel[sel != 1]} # any constrained model is preferred to the unconstrained fit.
      Fit <- Fit[sel]
      Covar <- Covar[sel]
      Covar.ok <- Covar.ok[sel]
      CrossIndex <- CrossIndex[sel]
      Lambda <- Lambda[sel]
    }

    win <- which.min(CrossIndex)
    fit <- Fit[[win]]; covar <- Covar[[win]]; covar.ok <- Covar.ok[win]; lambda <- Lambda[win]
  }


  fit$covar <- covar
  fit$covar.ok <- covar.ok
  fit$done <- done
  fit$warn <- FALSE
  fit$count <- count
  fit$lambda <- lambda
  fit
}

#' @export
qc.control <- function(maxTry = 25, trace = FALSE, lambda = NULL){
  if(!is.numeric(maxTry) || maxTry < 1 || !is.finite(maxTry))
    {stop("'maxTry' must be a positive integer, maxTry >= 1 is requested")}
  if(!is.logical(trace)){stop("'trace' must be TRUE/FALSE")}
  if(!is.null(lambda)){
    if(length(lambda) != 1){stop("'lambda' must be a scalar")}
    if(is.na(lambda) || lambda <= 0){stop("'lambda' must be a positive value")}
  }
  
  list(maxTry = maxTry, trace = trace, lambda = lambda)
}













