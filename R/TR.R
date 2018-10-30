TR <-
function(object, x){
	n <- dim(x)[1]
	p <- dim(x)[2]
	if (object$standardize){
		if (is.null(n)||(n<2)){
			stop("x must be a matrix with at least 2 observations!")
		}
	  if (n < p){
	    stop("Number of observations in x must be larger or at least the same as its dimension!")
	  }
	}
	if (class(object)=="MPL" ){
		## Standardize covariates
		if (object$standardize){
			x <- t(t(x)-colMeans(x))
			cholesky <- t(chol(cov(x)))
			x <- t(solve(cholesky, t(x)))
		}
		return (as.vector(((cbind(1, x) %*% object$beta.est) >0)+0))		
	}
	else {
		stop("Object must be fitted by MPL!")
	}
}
