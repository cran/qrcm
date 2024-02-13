# AUXILIARY FUNCTIONS FOR iqr


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# iqr auxiliary functions ###################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
summary.iqr <- function(object, p, cov = FALSE, ...){

	if(missing(p)){
		mf <- object$mf
		theta <- object$coefficients
		w <- attr(mf, "all.vars")$weights
  
		u <- sqrt(diag(object$covar))
		u <- matrix(u, q <- nrow(theta), r <- ncol(theta))
		dimnames(u) <- dimnames(theta)
		test <- (if(q*r == 1) NULL else iqr.waldtest(object))
		out <- list(converged = object$converged, n.it = object$n.it,
			coefficients = theta, se = u, 
			test.x = test$test.x, test.p = test$test.p)

		out$obj.function <- object$obj.function		
		out$n <- nrow(object$mf)
		out$free.par <- sum(object$s != 0)
		out$type <- attr(attr(object$mf, "type"), "fittype")
		out$n.events <- attr(object$mf, "n.events")
	}
	else{
		out <- list()
		for(i in 1:length(p)){out[[i]] <- extract.p(object, p[i], cov)}
		names(out) <- paste("p =", p)
		attr(out, "nacoef") <- which(apply(object$s,1, function(v){all(v == 0)}))
	}
	out$call <- object$call
	class(out) <- "summary.iqr"
	out	
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
print.summary.iqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if(!is.null(x$coef)){

		nacoef <- which(x$coef == 0)
		x$coef[nacoef] <- x$se[nacoef] <- NA

		cat("converged:", x$converged, "\n")
		cat("n. of iterations:", x$n.it, "\n")
		cat("n. of free parameters:", x$free.par, "\n")
		cat("n. of observations:", x$n, "\n")

		type <- x$type
		N <- x$n.events
		
		if(type %in% c("ctiqr", "ciqr")){
		  cat("N. of events: ", paste(deparse(round(N)), sep = " ", collapse = " "), "\n", sep = "")
		}
		else if(type == "iciqr"){
		  cat("* Events: ", paste(deparse(round(N[[2]])), sep = " ", collapse = " "), "\n", sep = "")
		  cat("* Left-censored: ", paste(deparse(round(N[[3]])), sep = " ", collapse = " "), "\n", sep = "")
		  cat("* Right-censored: ", paste(deparse(round(N[[1]])), sep = " ", collapse = " "), "\n", sep = "")
		  cat("* Interval-censored: ", paste(deparse(round(N[[4]])), sep = " ", collapse = " "), "\n", sep = "")
		}	

		cat("\n")
		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Coefficients:\n")
		print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("Standard errors:\n")
		print.default(format(x$se, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		if(!is.null(x$test.x)){
			cat("Wald test for x:\n")

			printCoefmat(x$test.x, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n\n")
		}

		if(!is.null(x$test.p)){
			cat("Wald test for b(p):\n")
			printCoefmat(x$test.p, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}

		cat("\n")
		if(is.null(attr(x$obj.function, "message"))){m <- "Minimized loss function:"}
		else{m <- "Loss (not the function being minimized):"}
		cat(m, x$obj.function)
		cat("\n\n")

	}

	else{
		nacoef <- attr(x, "nacoef")
		for(j in 1:(length(x) - 1)){
			cat(paste(names(x)[j], "\n"))
			cat("\n")
			cat("Coefficients:\n")
			coe <- x[[j]]$coef; coe[nacoef,] <- NA
			printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, 
				P.values = TRUE, has.Pvalue = TRUE)
			cat("\n")

			if(!is.null(x[[j]]$cov)){
				cat("Covar:\n")
				print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
			}
			cat("\n\n")
		}
	}

	invisible(x)
}



#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
plot.iqr <- function(x, conf.int = TRUE, polygon = TRUE, which = NULL, ask = TRUE, ...){

	plot.iqr.int <- function(p,u,j,conf.int,L){
		beta <- u[[j]]$beta
		if(is.null(L$ylim)){
			if(conf.int){y1 <- min(u[[j]]$low); y2 <- max(u[[j]]$up)}
			else{y1 <- min(beta); y2 <- max(beta)}
			L$ylim <- c(y1,y2)
		}
		plot(p, u[[j]]$beta, xlab = L$xlab, ylab = L$ylab, main = L$labels[j], 
		  type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col, axes = L$axes, 
		  frame.plot = L$frame.plot, cex.lab = L$cex.lab, cex.axis = L$cex.axis)
		if(conf.int){
		  if(polygon){
		    yy <- c(u[[j]]$low, tail(u[[j]]$up, 1), rev(u[[j]]$up), u[[j]]$low[1])
		    xx <- c(p, tail(p, 1), rev(p), p[1])
		    polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
		  }
		  else{
		    points(p, u[[j]]$low, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		    points(p, u[[j]]$up, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		  }
		}
	}

	L <- list(...)
	if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
	if(is.null(L$lwd)){L$lwd <- 2}
	if(is.null(L$col)){L$col <- "black"}
	if(is.null(L$xlab)){L$xlab <- "p"}
	if(is.null(L$ylab)){L$ylab <- "beta(p)"}
	if(is.null(L$cex.lab)){L$cex.lab <- 1}
	if(is.null(L$cex.axis)){L$cex.axis <- 1}
	if(is.null(L$axes)){L$axes <- TRUE}
	if(is.null(L$frame.plot)){L$frame.plot <- TRUE}
	L$labels <- rownames(x$coefficients)
	q <- length(L$labels)
	L$labels <- c(L$labels, "qqplot")

	p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
	u <- predict.iqr(x, p = p, type = "beta", se = conf.int)

	if(!is.null(which) | !ask){
		if(is.null(which)){which <- 1:q}
		for(j in which){plot.iqr.int(p,u,j,conf.int,L)}
	}
	else{
		pick <- 1
		while(pick > 0 && pick <= q + 1){
			pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
			if(pick > 0 && pick <= q){plot.iqr.int(p,u,pick,conf.int,L)}
			else if(pick == q + 1){
				KM <- attr(attr(x$mf, "CDFs"), "km")
				plot(KM$time, KM$cdf, pch = 20, cex = 0.5, 
				xlim = c(0,1), ylim = c(0,1), 
				ylab = "U(0,1) quantiles", xlab = "fitted CDF quantiles")
				if(!is.null(KM$low)){ # there is no CI if the data are IC
				  points(KM$time, KM$low, pch = ".")
				  points(KM$time, KM$up, pch = ".")
				}
				abline(0,1)
			}
		}
	}
}



#############################################################################################################
#############################################################################################################
#############################################################################################################



# predict function.
# p: default to percentiles for type = "beta". No default for "fitted". Ignored for "CDF".
# se: ignored for type = "CDF"
# x: only for type = "CDF" or type = "fitted"
# y: only for type = "CDF"

#' @export
predict.iqr <- function(object, type = c("beta", "CDF", "QF", "sim"), newdata, p, se = TRUE, ...){

	if(is.na(match(type <- type[1], c("beta", "CDF", "QF", "sim")))){stop("invalid 'type'")}
	if(type == "beta"){
		if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
		return(pred.beta(object, p, se))
	}

	mf <- object$mf
	mt <- terms(mf)
	miss <- attr(mf, "na.action")
	fittype <- attr(attr(mf, "type"), "fittype")
	nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	xlev <- .getXlevels(mt, mf)
	contr <- attr(mf, "contrasts")
	
	if(fittype == "iciqr" && missing(newdata) & type == "CDF")
	  {stop("with type = 'CDF' and interval-censored data, 'newdata' is always required \n and must include a variable named 'time' at which to compute predictions")}
	if(!missing(newdata)){

	  if(type == "CDF"){
	    yn <- as.character(if(fittype == "ctiqr") mt[[2]][[3]] 
              else if(fittype == "ciqr") mt[[2]][[2]] else if(fittype == "iqr") mt[[2]] else "time")
	    if(is.na(ind <- match(yn, colnames(newdata)))){
	      if(fittype == "icqr"){stop("with type = 'CDF' and interval-censored data, 'newdata' is always required \n and must include a variable named 'time' at which to compute predictions")}
	      else{stop("for 'type = CDF', 'newdata' must contain the y-variable")}
	    }
	    
	    # To have all the variables that WILL appear in all.vars(mt) but are not necessary for prediction
	    if(fittype == "ciqr"){newdata[,as.character(mt[[2]][[3]])] <- 1}
	    if(fittype == "ctiqr"){newdata[,as.character(mt[[2]][[4]])] <- 1
              newdata[,as.character(mt[[2]][[2]])] <- -Inf}
	    if(fittype == "iciqr"){mt <- update(mt, time ~ .)}
	  }
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}

		mf <- model.frame(mt, data = newdata, xlev = xlev)

		if(nrow(mf) == 0){
			nr <- nrow(newdata)
			if(type == "CDF"){
				out <- data.frame(matrix(NA,nr,3))
				colnames(out) <- c("log.f", "log.F", "log.S")
				rownames(out) <- rownames(newdata)
			}
			else if(type == "QF"){
				out <- data.frame(matrix(NA,nr,length(p)))
				colnames(out) <- paste("p",p, sep = "")
				rownames(out) <- rownames(newdata)
			}
			else{out <- rep.int(NA, nr)}
			return(out)
		}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
	}
	x <- model.matrix(mt, mf, contrasts.arg = contr)

	if(type == "CDF"){
		bfun <- attr(object$mf, "bfun")
		y <- cbind(model.response(mf))[,1 + (fittype == "ctiqr")]
		Fy <- p.bisec(object$coefficients, y,x, bfun)
		b1 <- apply_bfun(bfun, Fy, "b1fun")
		fy <- 1/c(rowSums((x%*%object$coefficients)*b1))
		# fy[attr(Fy, "out")] <- 0
		if(any(fy < 0)){warning("some PDF values are negative (quantile crossing)")}
		CDF <- PDF <- NULL
		CDF[nomiss] <- Fy
		PDF[nomiss] <- fy
		CDF[miss] <- PDF[miss] <- NA
		out <- data.frame(CDF = CDF, PDF = PDF)
		rownames(out)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(out)[miss] <- names(miss)}
		return(out)
	}

	else if(type == "QF"){
		if(missing(p)){stop("please indicate the value(s) of 'p' to compute x*beta(p)")}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}

		fit <- se.fit <- matrix(, length(c(miss,nomiss)), length(p))
		for(j in 1:length(p)){
			fit.beta <- extract.p(object,p[j], cov = se)
			fit[nomiss,j] <- x%*%cbind(fit.beta$coef[,1])
			if(se){se.fit[nomiss,j] <- sqrt(diag(x%*%fit.beta$cov%*%t(x)))}
		}
		fit <- data.frame(fit)
		colnames(fit) <- paste("p",p, sep = "")
		rownames(fit)[nomiss] <- rownames(mf)
		if(!is.null(miss)){rownames(fit)[miss] <- names(miss)}
		if(se){
			se.fit <- data.frame(se.fit)
			colnames(se.fit) <- paste("p",p, sep = "")
			rownames(se.fit)[nomiss] <- rownames(mf)
			if(!is.null(miss)){rownames(se.fit)[miss] <- names(miss)}
			return(list(fit = fit, se.fit = se.fit))
		}
		else{return(fit)}
	}	
	else{
		p <- runif(nrow(x))
		beta <- apply_bfun(attr(object$mf, "bfun"), p, "bfun")%*%t(object$coefficients)
		y <- NULL; y[nomiss] <- rowSums(beta*x); y[miss] <- NA
		return(y)
	}
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
terms.iqr <- function(x, ...){attr(x$mf, "terms")}
#' @export
model.matrix.iqr <- function(object, ...){
  mf <- object$mf
  mt <- terms(mf)
  model.matrix(mt, mf)
}
#' @export
vcov.iqr <- function(object, ...){object$covar}
#' @export
nobs.iqr <- function(object, ...){nrow(object$mf)}


#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
print.iqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

  
  #####
  
  cat("\n")
  cat("Number of free parameters:", sum(x$s), "\n")
  cat("Number of observations:", nrow(x$mf), "\n")

  type <- attr(attr(x$mf, "type"), "fittype")
  N <- attr(x$mf, "n.events")
  if(type %in% c("ctiqr", "ciqr")){
    cat("N. of events: ", paste(deparse(round(N)), sep = " ", collapse = " "), "\n", sep = "")
  }
  else if(type == "iciqr"){
    cat("* Events: ", paste(deparse(round(N[[2]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("* Left-censored: ", paste(deparse(round(N[[3]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("* Right-censored: ", paste(deparse(round(N[[1]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("* Interval-censored: ", paste(deparse(round(N[[4]])), sep = " ", collapse = " "), "\n", sep = "")
  }
  
  
  #####
  
  cat("\n")
	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	invisible(x)
}


