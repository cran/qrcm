


#' @export
test.fit.iqr <- function(object, R = 100, zcmodel = 1, trace = FALSE, ...){
  
  mf <- object$mf
  s <- object$s
  type <- attr(mf, "type")
  bfun <- attr(mf, "internal.bfun")
  bfun2 <- attr(mf, "bfun")
  theta <- attr(mf, "theta")
  theta2 <- object$coef
  Q0 <- attr(mf, "Q0"); chi0 <- qchisq(0.999, df = sum(s))

  statsy <- attr(object$mf, "stats")$y
  M <- 10/(statsy$M - statsy$m)
  V <- attr(mf, "all.vars"); U <- attr(mf, "all.vars.unscaled")
  x <- V$X; xw <- V$Xw; y <- V$y; z <- V$z0; d <- V$d; w <- V$weights; x2 <- U$X
  CDFs <- attr(mf, "CDFs")
  n <- N <- nrow(mf)
  q <- ncol(x)
  rho <- rep.int(1,n)
  maxit <- 10 + 2*sum(s)

  if(type == "ciqr"){
    Tc <- trans(z = -Inf, y = y, d = 1 - d, w = w, type = type)
    dat <- data.frame(z = -Inf, y = Tc$f(y), d = 1 - d, x = x)
    mc <- findagoodestimator(dat,w)
  }
  else if(type == "ctiqr"){
    minz <- min(z[z != -Inf])
    z <- pmax(z, minz)
    if(zcmodel == 1){u <- y - z; zz <- rep.int(0,n)}
    else if(zcmodel == 2){u <- y; zz <- z}
    else{stop("invalid 'zcmodel'")}


    Tz <- trans(z = -y, y = -z, d = rep.int(1,n), w = w, type = type)
    Tc <- trans(z = zz, y = u, d = 1 - d, w = w, type = type)

    dat.z <- data.frame(z = Tz$f(-y), y = Tz$f(-z), d = 1, x = x)
    dat.c <- data.frame(z = Tc$f(zz), y = Tc$f(u), d = 1 - d, x = x)

    # this should NOT happen. But sometimes it does, as the spline is not perfectly monotone
    dat.z$z[dat.z$z >= dat.z$y] <- -Inf 
    dat.c$z[dat.c$z >= dat.c$y] <- -Inf

    mz <- findagoodestimator(dat.z,w)
    mc <- findagoodestimator(dat.c,w)

    alphax <- alpha(object, mz,mc, zcmodel = zcmodel, Tc = Tc, Tz = Tz) # pr(Y > Z | x)
    N <- round(n/mean(alphax)) # n. of observations to sample in order to obtain n final obs.
    rho <- 1/alphax # inclusion probability of each obs.
  }

  exclude <- (if(type == "iqr") 0 else 0.02)
  test0 <- test.unif.ct(z = CDFs$CDF.z, y = CDFs$CDF.y, d = d, w = w, type = type, exclude)
  if(R == 0){
    out <- cbind(test0*c(1,sum(w)), NA)
    rownames(out) <- c("Kolmogorov-Smirnov", "Cramer-Von Mises")
    colnames(out) <- c("statistic", "p-value")
    return(out)
  }
  
  
  test <- matrix(,R,2)

  if(trace){pb <- txtProgressBar(min = 0, max = R, style = 3)}
  
  for(b in 1:R){
    
    if(trace){setTxtProgressBar(pb, b)}

    # x
    
    id <- sample.int(n, size = N, replace = TRUE, prob = rho)
    xb <- x[id,, drop = FALSE]
    xwb <- xw[id,, drop = FALSE]
    x2b <- x2[id,, drop = FALSE]
    wb <- w[id]
    
    # t
    beta <- tcrossprod(apply_bfun(bfun2, runif(N), "bfun"), theta2)
    tb <- yb <- (.rowSums(beta*x2b, N, q) - statsy$m)*M
    
    
    # z,c,y,d
    sim.pch <- getFromNamespace("sim.pch", ns = "pch")
    if(type == "ciqr"){
      cb <- (if(!is.null(mc)) sim.pch(mc, x = mc$x[id,,drop = FALSE], method = "s") else Inf)
      cb <- Tc$finv(cb)
      yb <- pmin(cb, tb)
      db <- (tb <= cb)
    }
    if(type == "ctiqr"){
      zb <- sim.pch(mz, x = mz$x[id,,drop = FALSE], method = "s")
      zb <- Tz$finv(zb)
      zb <- (minz - zb + abs(minz + zb))/2 # pmax(-zb, minz)
      cb <- (if(!is.null(mc)) sim.pch(mc, x = mc$x[id,,drop = FALSE], method = "s") else Inf)
      cb <- Tc$finv(cb)
      if(zcmodel == 1){cb <- cb + zb}
      yb <- pmin(cb, tb)
      db <- (tb <= cb)
      nb <- length(obs <- which(yb > zb))
      yb <- yb[obs]; zb <- zb[obs]; wb <- wb[obs]; db <- db[obs]
      xb <- xb[obs,, drop = FALSE]; xwb <- xwb[obs,, drop = FALSE]

      bfun$BB1 <- t(matrix(bfun$BB1[1,], ncol(bfun$BB1), nb))
    }
    
    # fit the model. Use small eps0: the starting points are generally less good than
    # those from start.theta!
   
    eps0 <- 0.001
    eeTol <- 1
    safeit <- 5
    chitol <- 1.5

    for(i in 1:5){
      
      fit.b <- iqr.newton(theta = theta, y = yb, z = zb, d = db, X = xb, Xw = xwb, 
                        bfun = bfun, s = s, type = type, 
                        maxit = maxit, tol = 1e-6, safeit = safeit, eps0 = eps0)
    
      thetadiff <- c(theta - fit.b$coef)[c(s == 1)]
      chi <- t(thetadiff)%*%Q0%*%cbind(thetadiff)
      
      check1 <- (fit.b$rank > 0)
      check2 <- (max(abs(fit.b$ee)) < eeTol)
      check3 <- (chi < chi0*chitol)
      
      if(fit.ok <- (check1 && check2 && check3)){break}
      eeTol <- eeTol * 2
      eps0 <- eps0/2
      safeit <- safeit + 2
      chitol <- chitol*2
    }

    if(fit.ok){
      test[b,] <- test.unif.ct(z = fit.b$pz, y = fit.b$py, 
          d = db, w = wb, type = type, exclude)
    }
  }
  if(trace){close(pb)}

  out <- cbind(test0*c(1,sum(w)), 
    c(
      mean(test[,1] >= test0[1], na.rm = TRUE), 
      mean(test[,2] >= test0[2], na.rm = TRUE)
    ))
  rownames(out) <- c("Kolmogorov-Smirnov", "Cramer-Von Mises")
  colnames(out) <- c("statistic", "p-value")

  out
}




#############################################################################################################
#############################################################################################################
#############################################################################################################


# kaplan-meier estimator for ct data. It returns time, cdf, n.event, cens, lower, upper.
# If exclude != NULL, it will exclude the indicated proportion of events (and of course all non-events
# in the same range). With Survival 2.41, a bug caused the function to stop with an error.
# To fix the bug, the function was split into "km" and "km.internal" as below. This was the
# only change from qrcm 2.0 to qrcm 2.1. In version 3.0, I added the option "timefix = FALSE" 
# which should remove the problem completely.

km <- function(z,y,d,w, type, exclude = NULL){

  km.internal <- function(z,y,d,w, type, exclude = NULL){
  
    if(type == "iqr"){m <- survfit(Surv(y, rep.int(1,length(y))) ~ 1, weights = w)}
    else if(type == "ciqr"){m <- survfit(Surv(y,d) ~ 1, weights = w)}
    else{
      m <- survfit(coxph(Surv(z,y,d) ~ 1, weights = pmax(w,1e-6), timefix = FALSE), type = "kaplan-meier")
    }

    out <- data.frame(time = m$time, cdf = 1 - m$surv, n.event = m$n.event, 
      low = 1 - m$upper, up = 1 - m$lower)
    out <- out[out$n.event > 1e-6,]

    if(type != "iqr" && !is.null(exclude)){
      u <- cumsum(out$n.event)
      u <- u/u[length(u)]
      u <- which(u >= 1 - exclude)
      out <- out[1:u[1],, drop = FALSE]
    }
    out
  }


  eps <- 1e-6
  n <- length(y)
  u <- (1:n)/n
  for(i in 1:10){
    out <- try(suppressWarnings(km.internal(z,y,d,w,type,exclude)), silent = TRUE)
    if(fit.ok <- !inherits(out, "try-error")){break}

    delta <- u*eps
    y <- y + delta
    z <- z - delta
    eps <- eps*2
  }
  out
}




#############################################################################################################
#############################################################################################################
#############################################################################################################



# ks and cvm statistic for uniformity with (possibly) censored and truncated data.
# Even if exclude = 0, the final censored observations are excluded anyway.
# Please avoid exclude = NULL, admitted but deprecated here.
test.unif.ct <- function(z,y,d,w, type, exclude = 0.05){
  
  y <- pmax(y,1e-8) # y = 0 causes problems, especially with ctiqr
  if(missing(w)){w <- NULL}
  KM <- suppressWarnings(km(z,y,d,w, type, exclude))

  n <- nrow(KM)
  hat.Fy <- KM$cdf
  Fy <- y <- KM$time # = punif(KM$time)
      
  # kolmogorov - smirnov
  
  DD <- Fy - hat.Fy
  ks <- max(abs(DD))
  
  # cramer - von mises
  
  Fy <- c(0, Fy)
  hat.Fy <- c(0, hat.Fy)
  y <- c(0,y)
  n <- n + 1
  
  U <- (hat.Fy - Fy)^2
  h <- y[2:n] - y[1:(n-1)]
  b1 <- U[1:(n-1)]
  b2 <- U[2:n]
  A <- (b1 + b2)*h/2
  cvm <- sum(A)

  ###
  
  c(ks = ks, cvm = cvm)	
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


# automatically finds a "good" estimator of a CDF. If there is no censoring (or very little censoring) on T,
# C may be (almost) completely censored. If this happens, I just return NULL.
findagoodestimator <- function(dat, w){
  if(sum(dat$d*w)/sum(w) < 0.05){return(NULL)}

  
  br <- max(5, min(15, round(nrow(dat)/30/(ncol(dat) - 3))))
  CDF <- suppressWarnings(pch::pchreg(
    Surv(z,y,d) ~ ., data = dat, weights = w, splinex = NULL, breaks = br)) 

  fit.ok <- (CDF$conv.status == 0)
  splx <- pch::splinex()
  br <- max(length(CDF$breaks),6)
  delta.br <- 1
  count <- 0
  while(!fit.ok){

    if(count == 32){break}
    br <- br + delta.br
    CDF <- suppressWarnings(pch::pchreg(
    Surv(z,y,d) ~ ., data = dat, weights = w, splinex = splx, breaks = br))
    count <- count + 1
  
    if(count == 5){br <- br - 5; delta.br <- -1}
    if(count == 10){br <- br + 5; delta.br <- 0; splx$df <- 3; splx$v <- 0.95}
    if(count == 11){delta.br <- 1}
    if(count == 16){br <- br - 5; delta.br <- -1}
    if(count == 21){br <- br + 5; delta.br <- 0; splx$df <- 1; splx$v <- 1}
    if(count == 22){delta.br <- 1}
    if(count == 27){br <- br - 5; delta.br <- -1}
    fit.ok <- (CDF$conv.status == 0)
  }
  # for directly applying pch:::sim.pch
  CDF$y <- attr(CDF$mf, "y")
  CDF$u <- attr(CDF$mf, "u")
  CDF$rangex <- attr(CDF$mf, "rangex")
  CDF$splinex <- attr(CDF$mf, "splinex")
  CDF
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


# computes P(Y > Z | x)
alpha <- function(obj, mz, mc, k = 98, zcmodel, Tc, Tz){
  
  p <- seq.int(0.02, 0.99, length.out = k) # avoid left tail for robustness
  t <- predict.iqr(obj, type = "QF", p = c(0.019, p), se = FALSE)
  t <- as.matrix(t); n <- nrow(t)
  staty <- attr(obj$mf, "stats")$y
  t <- (t - staty$m)/(staty$M - staty$m)*10
  
  Fz <- Sc <- matrix(,n,k+1)
  for(j in 1:(k + 1)){
    Fz[,j] <- quickpred(mz, y = Tz$f(-t[,j]), type = "SF")
    Sc[,j] <- (if(!is.null(mc) & zcmodel == 2) quickpred(mc, y = Tc$f(t[,j]), type = "SF") else 1)
  }
  Sc <- Sc[,-1]
  Fz0 <- Fz[,1]
  deltat <- t[,2:(k+1)] - t[,(1:k)]
  fz <- (Fz[,2:(k + 1)] - Fz[,1:k])/deltat
  
  r <- .rowSums(fz, n, k)
  r[r == 0] <- 1 # arbitrary, to avoid 0/0
  fz <- fz/r*(1 - Fz0)

  u <- cbind(Fz0, fz*Sc*t(matrix(1 - p, k,n)))
  u <- .rowSums(u, n, k + 1)
 
  pmax(0.1, pmin(u,1))
}



#############################################################################################################
#############################################################################################################
#############################################################################################################



quickpred <- function(obj, y, type = c("PDF", "SF")){

  type <- type[1]
  Breaks <- obj$breaks
  k <- attr(Breaks, "k")
  h <- attr(Breaks, "h")
  lambda <- obj$lambda
  Lambda <- obj$Lambda
  u <- attr(obj$mf, "u")
  end.y <- u(y); y <- (y + Breaks[1] - 1 + abs(y - Breaks[1] + 1))/2 # pmax(y, Breaks[1] - 1)
  n <- length(y)
  t <- y - c(0,Breaks)[end.y]
  ind <- cbind(1:n,end.y)
  lambda <- cbind(0,lambda)[ind]
  Lambda <- cbind(0,0,Lambda)[ind] + lambda*t
  SF <- exp(-Lambda)
  if(type == "SF"){return(SF)}
  else{return(lambda*SF)}
}



#############################################################################################################
#############################################################################################################
#############################################################################################################


# Transform a censored, truncated variable for a better prediction with pchreg.
# A common problem is that y may have few unique values, which may cause bad starting points
# in ctiqr.internal. To remove this problem, I add to y a small quantity that guarantees that
# length(unique(y)) = length(y).

trans <- function(z,y,d,w, type){

  if(all(d == 0)){return(list(f = I, finv = I))}
  
  # smoothing y
  yy <- sort(unique(y)); m <- length(yy)
  eps <- min(yy[2:m] - yy[1:(m - 1)])/10
  y <- y + (rank(y, ties.method = "first")/length(y)*eps)
  
  hatF <- km(z,y,d,w, type, exclude = NULL)
    
  n <- nrow(hatF)
  F0 <- hatF[1,]
  F1 <- hatF[n,]
  hatF <- hatF[2:(n - 1),]
  hatF <- hatF[!(hatF$cdf %in% 0:1),]
  
  n <- nrow(hatF)
  F0$cdf <- hatF$cdf[1]/10; F0$time <- F0$time - 1
  F1$cdf <- 1 - (1 - hatF$cdf[n])/10; F1$time <- F1$time + 1
  hatF <- rbind(F0, hatF, F1)
  
  n <- n + 2
  tstar <- -log(1 - hatF$cdf)
  t <- hatF$time
  tstar <- c(tstar[1] - 1, tstar, tstar[n] + 1)
  t <- c(t[1] - 10, t, t[n] + 10)
  
  # Don't use "hyman" (which is actually faster), as out-of-range predictions may be not monotone
  f <- splinefun(t, tstar, method = "monoH.FC")
  finv <- splinefun(tstar, t, method = "monoH.FC")

  list(f = f, finv = finv)
}
