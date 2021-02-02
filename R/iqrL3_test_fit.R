


# ks statistic for H0: u,v ~ independent U(0,1)
ks <- function(u,v,id,w1,w2, K = 25){

  v_long <- v; v <- v[!duplicated(id)]
  N <- length(u); n <- length(v)
  k1 <- pmin(K, floor(sqrt(N)))
  k2 <- pmin(K, floor(sqrt(n)))
  
  W2 <- sum(w2)
  p1 <- (1:k1)/k1
  p2 <- (1:k2)/k2
  F <- p2%*%t(p1)
  
  Fhat <- NULL
  for(i2 in 1:k2){
    Fv <- sum(w2*(v <= p2[i2]))/W2
    ww1 <- w1*(v_long <= p2[i2])
    W1 <- sum(ww1)
    for(i1 in 1:k1){
      Fu_v <- sum(ww1*(u <= p1[i1]))/W1
      Fhat <- c(Fhat, Fu_v*Fv)
    }
  }
  Fhat <- t(matrix(Fhat, k1,k2))

  list(F = data.frame(Fhat = c(Fhat), F = c(F)), 
       ks = max(abs(Fhat - F), na.rm = TRUE))
}







#' @export
test.fit.iqrL <- function(object, R = 100, trace = FALSE, ...){
  
  mf1 <- object$mf.theta
  mf2 <- object$mf.phi
  s.theta <- object$s.theta
  s.phi <- object$s.phi

  bfun1 <- attr(mf1, "internal.bfun")
  bfun1.bis <- attr(mf1, "bfun")
  theta <- attr(mf1, "theta")
  theta.bis <- object$theta

  bfun2 <- attr(mf2, "internal.bfun")
  bfun2.bis <- attr(mf2, "bfun")
  phi <- attr(mf2, "phi")
  phi.bis <- object$phi
 
  statsy1 <- attr(mf1, "stats")$y
  M1 <- 10/(statsy1$M - statsy1$m)
  statsy2 <- attr(mf2, "stats")$y
  M2 <- 10/(statsy2$M - statsy2$m)

  V1 <- attr(mf1, "all.vars"); U1 <- attr(mf1, "all.vars.unscaled")
  V2 <- attr(mf2, "all.vars"); U2 <- attr(mf2, "all.vars.unscaled")
  x <- V1$X; xw <- V1$Xw; w1 <- V1$w; x.bis <- U1$X
  z <- V2$X; zw <- V2$Xw; w2 <- V2$w; z.bis <- U2$X
  u <- object$fit$u
  v <- object$fit$v

  n <- nrow(mf1)
  n.id <- nrow(mf2)
  id <- as.numeric(as.factor(mf1$'(id)'))
  matchID <- match(mf1$'(id)',mf2$'(id)')
  q1 <- ncol(x); q2 <- ncol(z)

  maxit.theta <- 5 + 2*sum(s.theta)
  maxit.phi <- 5 + 2*sum(s.phi)
  maxit <- object$n.it*(1 + !object$converged) + maxit.theta + maxit.phi
  safeit.theta <- 2 + sum(s.theta)
  safeit.phi <- 2 + sum(s.phi)

  test0 <- ks(u,v,id,w1,w2, K = 25)$ks
  if(R == 0){
    out <- cbind(test0, NA)
    rownames(out) <- "Kolmogorov-Smirnov"
    colnames(out) <- c("statistic", "p-value")
    return(out)
  }
  
  test <- NULL
  if(trace){pb <- txtProgressBar(min = 0, max = R, style = 3)}

  for(b in 1:R){

    if(trace){setTxtProgressBar(pb, b)}

    # y - alpha
    beta <- tcrossprod(apply_bfun(bfun1.bis, runif(n), "bfun"), theta.bis)
    y_alpha <- (.rowSums(beta*x.bis, n, q1) - statsy1$m)*M1

    # alpha
    gamma <- tcrossprod(apply_bfun(bfun2.bis, runif(n.id), "bfun"), phi.bis)
    alpha <- (.rowSums(gamma*z.bis, n.id, q2) - statsy2$m)*M2


    # fit the model
    eps0 <- 0.025 # more prudent than iqrL.iternal: "true" starting points may be far from the solution
    gTol <- 1
    for(i in 1:5){
      fit <- iqrL.fit(theta,phi, y_alpha + alpha[matchID], alpha, x,xw,z,zw, id, w1,w2, bfun1,bfun2, s.theta,s.phi, 
			maxit.theta + i, safeit.theta + i, maxit.phi + i, safeit.phi + i, eps0, tol = 1e-5, maxit)
      fit.ok <- (fit$fullrank & max(abs(fit$g)/sqrt(n)) < gTol)
      if(fit.ok){break}
      gTol <- gTol * 2
      eps0 <- eps0/2
    }
    if(fit.ok){test[b] <- ks(fit$A$u,fit$A$v[matchID],id,w1,w2, K = 25)$ks}
  }
  if(trace){close(pb)}
  out <- cbind(test0, mean(test >= test0, na.rm = TRUE))
  rownames(out) <- "Kolmogorov-Smirnov"
  colnames(out) <- c("statistic", "p-value")

  out
}

