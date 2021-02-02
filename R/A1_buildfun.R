# FUNCTIONS FOR MODEL BUILDING 

#############################################################################################################
#############################################################################################################
#############################################################################################################



#' @export
plf <- function(p, knots){ # basis of piecewise linear function
	if(is.null(knots)){return(cbind(b1 = p))}
	k <- length(knots <- sort(knots))
	ind1 <- 1
	ind2 <- NULL
	for(j in knots){ind1 <- cbind(ind1, (p > j))}
	for(j in k:1){ind2 <- cbind(ind2, ind1[,j] - ind1[,j+1])}
	ind2 <- cbind(ind1[,k+1], ind2)[,(k + 1):1, drop = FALSE]
	ind1 <- cbind(ind1,0)

	knots <- c(0,knots,1)
	a <- NULL
	for(j in 1:(k + 1)){
		a <- cbind(a, (p - knots[j])*ind2[,j] + (knots[j + 1] - knots[j])*ind1[,j + 1])
	}
	colnames(a) <- paste("b", 1:(k + 1), sep = "")
	attr(a, "knots") <- knots[2:(k+1)]
	a
}


#############################################################################################################
#############################################################################################################
#############################################################################################################




slp.basis <- function(k, intercept){ # shifted Legendre polynomials basis

	K <- k + 1

	# matrix a such that P%*%a is an orthogonal polynomial, P = (1, p, p^2, p^3, ...)

	a <- matrix(0, K,K)
	for(i1 in 0:k){
		for(i2 in 0:i1){
			a[i2 + 1, i1 + 1] <- choose(i1,i2)*choose(i1 + i2, i2)
		}
	}
	a[,seq.int(2, K, 2)] <- -a[,seq.int(2, K, 2)]
	a[seq.int(2, K, 2),] <- -a[seq.int(2, K, 2),]

	# a1 = first derivatives to be applied to P' = (0,1, p, p^2, ...)
	# A = first integral to be applied to PP = (p, p^2, p^3, p^4, ...)
	# AA = second integral to be applied to PPP = (p^2, p^3, p^4, p^5, ...)

	a1 <- A <- AA <- matrix(,K,K)
	for(j in 0:k){
		a1[j + 1,] <- a[j + 1,]*j
		A[j + 1,] <- a[j + 1,]/(j + 1)
		AA[j + 1,] <- A[j + 1,]/(j + 2)
	}

	if(!intercept){a[1,-1] <- A[1,-1] <- AA[1, -1] <- 0}
	out <- list(a = a, a1 = a1, A = A, AA = AA)
	attr(out, "k") <- k
	class(out) <- "slp.basis"
	out
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


#' @export
slp <- function(p, k = 3, intercept = FALSE){
	if((k <- round(k)) < 1){stop("k >= 1 is required")}
	P <- cbind(1, outer(p, seq_len(k), "^"))
	B <- P%*%slp.basis(k, intercept)$a
	colnames(B) <- paste("slp", 0:k, sep = "")
	B <- B[,-1, drop = FALSE]
	attr(B, "k") <- k
	class(B) <- "slp"
	B
}


#############################################################################################################
#############################################################################################################
#############################################################################################################


is.slp <- function(f){
	test.p <- seq(0,1,0.1)
	B <- model.matrix(f, data = data.frame(p = test.p, u = test.p, v = test.p))
	if(nrow(B) == 0){return(FALSE)}
	a <- attr(B, "assign")
	if(any(a > 1)){return(FALSE)}
	B <- B[,a == 1, drop = FALSE]
	k <- ncol(B)
	intercept <- FALSE
	if(any(B != slp(test.p, k = k, intercept))){
		intercept <- TRUE
		if(any(B != slp(test.p, k = k, intercept))){
			return(FALSE)
		}
	}
	out <- TRUE
	attr(out, "k") <- k
	attr(out, "intercept") <- intercept
	attr(out, "intB") <- any(a == 0)
	out
}

