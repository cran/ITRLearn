maximin <- function(B, c0)
# B is the matrix containing (estimated) coefficients in groupwise optimal ITR
# c0 is the estimated marginal treatment effect
{
	p <- dim(B)[1]
	G0 <- dim(B)[2]
	Cmat <- rep(0, G0)
	Amat <- matrix(1, 1, G0)
	Dmat <- crossprod(B, B)
	bvec <- 1
	rvec <- 0
	lvec <- rep(0, G0)
	uvec <- rep(1, G0)
	qp <- ipop(c=Cmat, H=Dmat, A=Amat, b=bvec, l=lvec, u=uvec, r=0)
	x0 <- primal(qp)
	t0 <- crossprod(x0, Dmat%*%x0)
	theta0 <- rep(0, p+1)
	if (t0<=0){
		theta0[1] <- (c0>0)+0
		I0 <- 1:G0
	}
	else{
		theta0[2:(p+1)] <- B %*% x0
		theta0[2:(p+1)] <- theta0[2:(p+1)]/sqrt(sum(theta0[2:(p+1)]^2))
		t0 <- min(crossprod(B, theta0[2:(p+1)]))
		theta0[2:(p+1)] <- theta0[2:(p+1)]*t0/sqrt(t0^2+c0^2)
		theta0[1] <- c0/sqrt(t0^2+c0^2)
		I0 <- which(abs(crossprod(B, theta0[2:(p+1)])-t0)<1e-16)
	}
	
	result <- list(beta.est=theta0, I0=I0)
	class(result) <- "MPL"
	return(result)
}