# AUXILIARY FUNCTIONS FOR iqrL


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
# iqrL auxiliary functions ##################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################




#' @export
summary.iqrL <- function(object, p, level, cov = FALSE, ...){

	o <- object
	o1 <- list(mf = o$mf.theta, coefficients = o$theta, s = o$s.theta, covar = o$covar.theta)
	o2 <- list(mf = o$mf.phi, coefficients = o$phi, s = o$s.phi, covar = o$covar.phi)
	if(missing(p)){

		theta <- o$theta
		v1 <- sqrt(diag(o$covar.theta))
		v1 <- matrix(v1, q1 <- nrow(theta), k1 <- ncol(theta))
		dimnames(v1) <- dimnames(theta)

		phi <- o$phi
		v2 <- sqrt(diag(o$covar.phi))
		v2 <- matrix(v2, q2 <- nrow(phi), k2 <- ncol(phi))
		dimnames(v2) <- dimnames(phi)

		test1 <- (if(q1*k1 == 1) NULL else iqr.waldtest(o1))
		test2 <- (if(q2*k2 == 1) NULL else iqr.waldtest(o2))

		out <- list(converged = o$converged, n.it = o$n.it,
			theta = theta, se.theta = v1, 
			phi = phi, se.phi = v2, 
			test.row.theta = test1$test.x, test.col.theta = test1$test.p,
			test.row.phi = test2$test.x, test.col.phi = test2$test.p
		)

		out$obj.function <- o$obj.function[1]
		out$n <- nrow(o$mf.theta)
		out$n.id <- nrow(o$mf.phi)
		out$free.par <- sum(o$s.theta) + sum(o$s.phi)
	}
	else{
		if(!is.atomic(level) || length(level) != 1 || !(level %in% 1:2)){
			stop("set level = 1 to summarize beta(u), and level = 2 to summarize gamma(v)")
		}
		oo <- (if(level == 1) o1 else o2)
		s <- (if(level == 1) o$s.theta else o$s.phi)
		out <- list()
		for(i in 1:length(p)){out[[i]] <- extract.p(oo, p[i], cov)}
		names(out) <- (if(level == 1) paste("u =", p) else paste("v =", p))
		attr(out, "nacoef") <- which(apply(s,1, function(h){all(h == 0)}))
		attr(out, "p") <- p
		attr(out, "is.cov") <- cov
		attr(out, "level") <- level
	}
	out$call <- o$call
	class(out) <- "summary.iqrL"
	out	
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
print.summary.iqrL <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	cat("######################", "\n")
	cat("######################", "\n\n")
	if(!is.null(x$converged)){

		natheta <- which(x$theta == 0)
		x$theta[natheta] <- x$se.theta[natheta] <- NA
		naphi <- which(x$phi == 0)
		x$phi[naphi] <- x$se.phi[naphi] <- NA

		cat("converged:", x$converged, "\n")
		cat("n. of iterations:", x$n.it, "\n")
		cat("n. of observations:", x$n, "\n")
		cat("n. of clusters:", x$n.id, "\n")
		cat("n. of free parameters (excl. alpha):", x$free.par, "\n\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Coefficients: theta\n")
		print.default(format(x$theta, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("Standard errors: theta\n")
		print.default(format(x$se.theta, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Coefficients: phi\n")
		print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")

		cat("Standard errors: phi\n")
		print.default(format(x$se.phi, digits = digits), print.gap = 2L, quote = FALSE)
		cat("\n")
		
		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Wald test for rows of theta:\n")
		if(!is.null(x$test.row.theta)){
			printCoefmat(x$test.row.theta, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}
		else{cat("(omitted - theta has a single row)")}
		cat("\n\n")

		cat("Wald test for columns of theta:\n")
		if(!is.null(x$test.col.theta)){
			printCoefmat(x$test.col.theta, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}
		else{cat("(omitted - theta has a single column)")}
		cat("\n\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Wald test for rows of phi:\n")
		if(!is.null(x$test.row.phi)){
			printCoefmat(x$test.row.phi, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}
		else{cat("(omitted - phi has a single row)")}
		cat("\n\n")

		cat("Wald test for columns of phi:\n")
		if(!is.null(x$test.col.phi)){
			printCoefmat(x$test.col.phi, digits = digits, signif.stars = TRUE, 
				signif.legend = FALSE, zap.ind = 2, tst.ind = 1, 
				P.values = TRUE, has.Pvalue = TRUE)
		}
		else{cat("(omitted - phi has a single column)")}
		cat("\n\n")

		cat("######################", "\n")
		cat("######################", "\n\n")

		cat("Minimized loss function:", x$obj.function)
		cat("\n")
		cat("\n")
	}

	else{
		level <- attr(x, "level")
		nacoef <- attr(x, "nacoef")
		is.cov <- attr(x, "is.cov")
		p <- attr(x, "p")
		coefname <- (if(level == 1) "beta" else "gamma")
		for(j in 1:(length(x) - 1)){

			cat("Coefficients:", paste0(coefname,"(",p[j],")"), "\n")
			coe <- x[[j]]$coef; coe[nacoef,] <- NA
			printCoefmat(coe, digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, 
				P.values = TRUE, has.Pvalue = TRUE)

			cat("\n")

			if(is.cov){
				cat("######################", "\n\n")
				cat("Covar:", paste0(coefname, "(",p[j],")"), "\n")
				print.default(format(x[[j]]$cov, digits = digits), print.gap = 2L, quote = FALSE)
			}
			cat("\n")
			cat("######################", "\n")
			cat("######################", "\n\n")
		}
	}

	invisible(x)
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
plot.iqrL <- function(x, conf.int = TRUE, polygon = TRUE, which = NULL, ask = TRUE, ...){

	plot.iqrL.int <- function(p,pred,j,conf.int,L, level){
  
		coe <- pred[[j]][,2]; low <- pred[[j]]$low; up <- pred[[j]]$up
		if(is.null(L$ylim)){
			if(conf.int){y1 <- min(low); y2 <- max(up)}
			else{y1 <- min(coe); y2 <- max(coe)}
			L$ylim <- c(y1,y2)
		}
		if(is.null(L$xlab)){L$xlab <- (if(level == 1) "u" else "v")}
		if(is.null(L$ylab)){L$ylab <- (if(level == 1) "beta(u)" else "gamma(v)")}
		L$labels <- (if(level == 1) L$labels1 else L$labels2)

		plot(p, coe, xlab = L$xlab, ylab = L$ylab, main = L$labels[j], 
		type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col)
		if(conf.int){

		  if(polygon){
		    yy <- c(low, tail(up, 1), rev(up), low[1])
		    xx <- c(p, tail(p, 1), rev(p), p[1])
		    polygon(xx, yy, col = adjustcolor(L$col, alpha.f = 0.25), border = NA)
		  }
		  else{
		    points(p, low, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		    points(p, up, lty = 2, lwd = L$lwd, type = "l", col = L$col)
		  }
		}
	}

	q1 <- nrow(x$theta)
	q2 <- nrow(x$phi)
	q <- q1 + q2
	level <- c(rep(1,q1), rep(2,q2))
	L <- list(...)
	if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
	if(is.null(L$lwd)){L$lwd <- 2}
	if(is.null(L$col)){L$col <- "black"}
	L$labels1 <- paste0("beta:", rownames(x$theta))
	L$labels2 <- paste0("gamma:", rownames(x$phi))
	L$labels <- c(L$labels1, L$labels2, "qqplot(u)", "qqplot(v)", "ppplot(u,v)")

	p <- seq.int(max(0.001,L$xlim[1]), min(0.999,L$xlim[2]), length.out = 100)
	pred1 <- predict.iqrL(x, level = 1, p = p, type = "coef", se = conf.int)
	pred2 <- predict.iqrL(x, level = 2, p = p, type = "coef", se = conf.int)

	if(!is.null(which) | !ask){
		if(is.null(which)){which <- 1:q}
		for(pick in which){
			pred <- (if(pick <= q1) pred1 else pred2)
			j <- (if(pick <= q1) pick else pick - q1)
			plot.iqrL.int(p,pred,j,conf.int,L,level[pick])
		}
	}
	else{
		pick <- 1
		while(pick > 0 && pick <= q + 3){
			pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
			if(pick > 0 && pick <= q){
				pred <- (if(pick <= q1) pred1 else pred2)
				j <- (if(pick <= q1) pick else pick - q1)
				plot.iqrL.int(p,pred,j,conf.int,L, level[pick])
			}
			else{
				u <- x$fit$u
				vlong <- x$fit$v
				v <- vlong[!duplicated(x$fit$id)]
				id <- x$fit$id
				w1 <- attr(x$mf.theta, "all.vars")$w
				w2 <- attr(x$mf.phi, "all.vars")$w
				f <- function(x){1 - x}
				if(pick == q + 1){
					KM <- survfit(Surv(u) ~ 1, weights = w1)
					plot(KM, xlim = c(0,1), ylab = "U(0,1) quantiles", xlab = "fitted u", fun = f)
					abline(0,1)
				}
				else if(pick == q + 2){
					KM <- survfit(Surv(v) ~ 1, weights = w2)
					plot(KM, xlim = c(0,1), ylab = "U(0,1) quantiles", xlab = "fitted v", fun = f)
					abline(0,1)
				}
				else if(pick == q + 3){
				  Fuv <- ks(u,vlong,id,w1,w2, K = 10)$F
					plot(Fuv$Fhat, Fuv$F, xlim = c(0,1), ylim = c(0,1), 
					     xlab = expression(hat(F)(u,v)), ylab = "F(u,v)")
					abline(0,1)
				}
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
predict.iqrL <- function(object, level, type = c("coef", "CDF", "QF", "sim"), newdata, p, se = FALSE, ...){

	if(is.na(match(type <- type[1], c("coef", "CDF", "QF", "sim")))){stop("invalid 'type'")}
	if(!is.atomic(level) || length(level) != 1 || !(level %in% 1:2)){
		if(type == "coef"){stop("set level = 1 to predict beta, and level = 2 to predict gamma")}
		if(type == "CDF"){stop("set level = 1 to predict U, and level = 2 to predict V")}
		if(type == "QF"){stop("set level = 1 to predict quantiles of y - alpha, and level = 2 to predict quantiles of alpha")}
		if(type == "sim"){stop("set level = 1 to simulate y - alpha, and level = 2 to simulate alpha")}
	}
	o <- object
	if(level == 1){oo <- list(mf = o$mf.theta, coefficients = o$theta, covar = o$covar.theta, y = o$fit$y_alpha)}
	else{oo <- list(mf = o$mf.phi, coefficients = o$phi, covar = o$covar.phi, y = o$fit$alpha[!duplicated(o$fit$id)])}
	predict_iqrL.internal(oo, level = level, type = type, newdata = newdata, p = p, se = se)
}

# slightly different from qrcm::predict.iqr
predict_iqrL.internal <- function(object, level, type = c("coef", "CDF", "QF", "sim"), newdata, p, se = FALSE, ...){

	pname <- (if(level == 1) "u" else "v")
	coename <- (if(level == 1) "beta" else "gamma")
	if(type == "coef"){
		if(missing(p)){p <- seq.int(0.01,0.99,0.01)}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
		out <- pred.beta(object, p, se)
		for(j in 1:length(out)){
		  names(out[[j]])[1:2] <- c(pname, coename)
		}
		return(out)
	}

	mf <- object$mf
	mt <- terms(mf)
	xlev <- .getXlevels(mt, mf)
	contr <- attr(mf, "contrasts")

	if(datain <- !missing(newdata)){

		if(type == "CDF"){
			yn <- (if(level == 1) "y_alpha" else "alpha")
			if(is.na(ind <- match(yn, colnames(newdata))))
			{stop("for 'type = CDF', 'newdata' must contain the response variable 
				('y_alpha', if level = 1, and 'alpha', if level = 2)")}
			y <- newdata[,ind]
		}
		mt <- delete.response(mt)
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}

		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){stop("no non-missing values in mf")}
		if(type == "CDF" && any(miss <- attr(mf, "na.action"))){y <- y[-miss]}
	}
	else{if(type == "CDF"){y <- object$y}}

	x <- model.matrix(mt, mf, contrasts.arg = contr)

	if(type == "CDF"){
		bfun <- attr(object$mf, "bfun")
		Fy <- p.bisec(object$coefficients, y,x, bfun)
		b1 <- apply_bfun(bfun, Fy, "b1fun")
		fy <- 1/c(rowSums((x%*%object$coefficients)*b1))
		# fy[attr(Fy, "out")] <- 0
		if(any(fy < 0)){warning("some PDF values are negative (quantile crossing)")}
		out <- data.frame(CDF = Fy, PDF = fy)
		rownames(out) <- (if(!datain & level == 2) mf[,"(id)"] else rownames(mf))
		return(out)
	}

	else if(type == "QF"){
		if(missing(p)){stop("please indicate the value(s) of 'p' to compute quantiles")}
		if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}

		fit <- se.fit <- matrix(, nrow(mf), length(p))
		colnames(fit) <- colnames(se.fit) <- paste0(pname,p)
		rownames(fit) <- rownames(se.fit) <- (if(!datain & level == 2) mf[,"(id)"] else rownames(mf))
		for(j in 1:length(p)){
			fit.beta <- extract.p(object,p[j], cov = se)
			fit[,j] <- x%*%cbind(fit.beta$coef[,1])
			if(se){se.fit[,j] <- sqrt(diag(x%*%fit.beta$cov%*%t(x)))}
		}
		fit <- data.frame(fit)
		if(se){
			se.fit <- data.frame(se.fit)
			return(list(fit = fit, se.fit = se.fit))
		}
		else{return(fit)}
	}	
	else{
		p <- runif(nrow(x))
		beta <- apply_bfun(attr(object$mf, "bfun"), p, "bfun")%*%t(object$coefficients)
		y <- rowSums(beta*x)
		names(y) <- (if(!datain & level == 2) mf[,"(id)"] else rownames(mf))
		return(y)
	}
}


#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
terms.iqrL <- function(x, ...){
	list(
		terms1 = attr(x$mf.theta, "terms"),
		terms2 = attr(x$mf.phi, "terms")
	)
}
#' @export
model.matrix.iqrL <- function(object, ...){
  mf1 <- object$mf.theta
  mt1 <- terms(mf1)
  x1 <- model.matrix(mt1, mf1)

  mf2 <- object$mf.phi
  mt2 <- terms(mf2)
  x2 <- model.matrix(mt2, mf2)

  list(x1 = x1, x2 = x2)
}
#' @export
vcov.iqrL <- function(object, ...){list(covar.theta = object$covar.theta, covar.phi = object$covar.phi)}

#' @export
nobs.iqrL <- function(object, ...){list(n = nrow(object$mf.theta), n.id = nrow(object$mf.phi))}



#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
print.iqrL <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

	cat("Coefficients: theta\n")
	print.default(format(x$theta, digits = digits), print.gap = 2L, quote = FALSE)
	
	cat("\n")
	cat("Coefficients: phi\n")
	print.default(format(x$phi, digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	cat("Minimized Loss function:\n")
	print.default(format(as.numeric(x$obj.function[1]), digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	invisible(x)
}

