MPL.fit <-
function(y, x.tau, a, g=NULL, x.h=NULL, x.pi=NULL, method=c("Q", "A"),
			bootstrap=FALSE, random=FALSE, control=MPL.control())
{
	# a is binary or not
	if (any((a!=0) & (a!=1)))
		stop("Treatment must be binary variable!")
	# number of observations	
	if (length(y) < 1)
		stop("empty model")
	# A learning is used by default
	if (missing(method))
		method <- "A"	
	# bootstrap cannot be used if the baseline or propensity score is prespecified
	if (bootstrap&&((!is.null(control$h.est))||(!is.null(control$pi.est))))
		stop("Bootstrap cannot be used if estimators for the baseline or the propensity score is prespecified!")
	# x.tau
	x.tau <- as.matrix(x.tau)
	# x.h
	if (is.null(x.h)){
		x.h <- x.tau
	}
	else{
		x.h <- as.matrix(x.h)
	}
	# x.pi
	if (is.null(x.pi)){
		x.pi <- x.tau
	}
	else{
		x.pi <- as.matrix(x.pi)
	}
	
	# Q learning	
	if (method=="Q"){
		result <- Q.fit(y=y, x.tau=x.tau, a=a, g=g, x.h=x.h, control=control)
	}	
	# A learning
	else{
		result <- A.fit(y=y, x.tau=x.tau, a=a, g=g, x.h=x.h, x.pi=x.pi, random=random, control=control)
	}
	
	# maximin projection learning
	if (is.null(result$B)||is.null(result$c0)){
		result$beta.est <- NULL
	}
	else{
		result$beta.est <- maximin(result$B, result$c0)$beta.est
	}
	
	# random or not
	pi.est <- control$pi.est
	h.est <- control$h.est
	boot.sample <- control$boot.sample
	if (bootstrap){
		if (is.null(g)||(length(unique(g))==1)){
			Theta.tau.boot <- matrix(0, 1+dim(x.tau)[2], boot.sample)
			if (is.null(h.est)){
				Theta.h.boot <- matrix(0, 1+dim(x.h)[2], boot.sample)
			}
			if (is.null(pi.est)&&(!random)&&(method=="A")){
				Theta.pi.boot <- matrix(0, 1+dim(x.pi)[2], boot.sample)
			}
			if (method=="Q"){
				for (b in 1:boot.sample){
				  Index <- sample(1:length(y), length(y), replace=TRUE)
					y.boot <- y[Index]
					a.boot <- a[Index]
					x.tau.boot <- x.tau[Index,]
					x.h.boot <- x.h[Index,]
					md <- Q.fit(y=y.boot, x.tau=x.tau.boot, a=a.boot, g=g, x.h=x.h.boot, control=control)
					Theta.tau.boot[,,b] <- md$Theta.tau.est
					if (is.null(h.est)){
						Theta.h.boot[,,b] <- md$Theta.h.est
					}
				}
			}
			else{
				for (b in 1:boot.sample){
				  Index <- sample(1:length(y), length(y), replace=TRUE)
					y.boot <- y[Index]
					a.boot <- a[Index]
					x.tau.boot <- x.tau[Index,]
					x.h.boot <- x.h[Index,]
					x.pi.boot <- x.pi[Index,]
					md <- A.fit(y=y.boot, x.tau=x.tau.boot, x.pi=x.pi.boot, a=a.boot, g=g, x.h=x.h.boot, random=random, control=control)
					Theta.tau.boot[,,b] <- md$Theta.tau.est
					if (is.null(h.est)){
						Theta.h.boot[,,b] <- md$Theta.h.est
					}
					if (is.null(pi.est)&&(!random)){
						Theta.pi.boot[,,b] <- md$Theta.pi.est
					}
				}
			}
			result$Theta.tau.boot <- Theta.tau.boot
			if (is.null(h.est)){
			  result$Theta.h.boot <- Theta.h.boot
			}
			if (is.null(pi.est)&&(!random)){
			  result$Theta.pi.boot <- Theta.pi.boot
			}
		}
		else{
			G <- length(unique(g))
			Theta.tau.boot <- array(0, c(dim(x.tau)[2]+1, G, boot.sample))
			beta.est.boot <- matrix(0, dim(x.tau)[2]+1, boot.sample)
			if (is.null(h.est)){
				Theta.h.boot <- array(0, c(dim(x.h)[2]+1, G, boot.sample))
			}
			if (is.null(pi.est)&&(!random)&&(method=="A")){
				Theta.pi.boot <- array(0, c(dim(x.pi)[2]+1, G, boot.sample))
			}
			if (method=="Q"){
				for (b in 1:boot.sample){
				  Index <- rep(0, length(y))
				  for (j in 1:G){
				    Index[g==unique(g)[j]] <- sample(which(g==unique(g)[j]), sum(g==unique(g)[j]), replace=TRUE)
				  }
					y.boot <- y[Index]
					a.boot <- a[Index]
					x.tau.boot <- x.tau[Index,]
					x.h.boot <- x.h[Index,]
					md <- Q.fit(y=y.boot, x.tau=x.tau.boot, a=a.boot, g=g, x.h=x.h.boot, control=control)
					Theta.tau.boot[,,b] <- md$Theta.tau.est
					if (is.null(h.est)){
						Theta.h.boot[,,b] <- md$Theta.h.est
					}
					beta.est.boot[,b] <- maximin(md$B, md$c0)$beta.est
				}
			}
			else{
				for (b in 1:boot.sample){
				  Index <- rep(0, length(y))
				  for (j in 1:G){
				    Index[g==unique(g)[j]] <- sample(which(g==unique(g)[j]), sum(g==unique(g)[j]), replace=TRUE)
				  }
					y.boot <- y[Index]
					a.boot <- a[Index]
					x.tau.boot <- x.tau[Index,]
					x.h.boot <- x.h[Index,]
					x.pi.boot <- x.pi[Index,]
					md <- A.fit(y=y.boot, x.tau=x.tau.boot, x.pi=x.pi.boot, a=a.boot, g=g, x.h=x.h.boot, random=random, control=control)
					Theta.tau.boot[,,b] <- md$Theta.tau.est
					if (is.null(h.est)){
						Theta.h.boot[,,b] <- md$Theta.h.est
					}
					if (is.null(pi.est)&&(!random)){
						Theta.pi.boot[,,b] <- md$Theta.pi.est
					}
					beta.est.boot[,b] <- maximin(md$B, md$c0)$beta.est
				}
			}
			result$Theta.tau.boot <- Theta.tau.boot
			result$beta.boot <- beta.est.boot
			if (is.null(h.est)){
			  result$Theta.h.boot <- Theta.h.boot
			}
			if (is.null(pi.est)&&(!random)){
			  result$Theta.pi.boot <- Theta.pi.boot
			}
		}
	}
	
	class(result) <- "MPL"
	if (is.null(g)||(length(unique(g))==1)){
		result$standardize <- FALSE
	}
	else{
		result$standardize <- TRUE
	}	
	return(result)
}

Q.fit <-
function(y, x.tau, a, g=NULL, x.h=NULL, control=MPL.control())
{
	pi.est <- control$pi.est
	h.est <- control$h.est
	boot.sample <- control$boot.sample
	if (!is.null(h.est)){
		Theta.h.est <- NULL
		y <- y-h.est
		if (is.null(g)||(length(unique(g))==1)){
			x0 <- a*cbind(1,x.tau)
			Theta.est <- as.vector(lm(y~x0-1))
			B <- NULL
			c0 <- NULL
		}
		else{
			p <- dim(x.tau)[2]
			G <- length(unique(g))
			Theta.est <- matrix(0, p+1, G)
			mu <- matrix(0, p, G)
			cholesky <- array(0, c(p, p, G))
			x0 <- matrix(0, length(y), p*G+1)
			x0[,1] <- a
			for (j in 1:G){
				mu[,j] <- colMeans(x.tau[g==unique(g)[j],])
				cholesky[,,j] <- t(chol(cov(x.tau[g==unique(g)[j],])))
				x.tau[g==unique(g)[j],] <- t(t(x.tau[g==unique(g)[j],])-mu[,j])
				x.tau[g==unique(g)[j],] <- t(solve(cholesky[,,j], t(x.tau[g==unique(g)[j],])))
			  x0[,((j-1)*p+2):(j*p+1)] <- x.tau * ((g==unique(g)[j])+0) * a
			}
			coef0 <- as.vector(coef(lm(y~x0-1)))
			Theta.est[1,] <- coef0[1]
			c0 <- coef0[1]
			Theta.est[2:(p+1),] <- matrix(coef0[-1], p, G)
			B <- matrix(coef0[-1], p, G)
			for (j in 1:G){
				Theta.est[-1,j] <- solve(t(cholesky[,,j]), Theta.est[-1,j])
				Theta.est[1,j] <- Theta.est[1,j]-crossprod(Theta.est[-1,j], mu[,j])
			}
		}
	}
	
	else{
		p <- dim(x.tau)[2]
		if (is.null(g)||(length(unique(g))==1)){
			x0 <- cbind(a*cbind(1,x.tau), 1, x.h)
			Theta.est <- as.vector(lm(y~x0-1))
			Theta.h.est <- Theta.est[-(1:(p+1))]
			Theta.est <- Theta.est[1:(p+1)]
			h.est <- cbind(1,x.h)%*%Theta.h.est
			B <- NULL
			c0 <- NULL
		}
		else{
			G <- length(unique(g))
			p.h <- dim(x.h)[2]
			Theta.est <- matrix(0, p+1, G)
			Theta.h.est <- matrix(0, p.h+1, G)
			x0 <- matrix(0, length(y), p*G+G*(p.h+1)+1)
			x0[,1] <- a
			mu <- matrix(0, p, G)
			cholesky <- array(0, c(p, p, G))
			for (j in 1:G){
				mu[,j] <- colMeans(x.tau[g==unique(g)[j],])
				cholesky[,,j] <- t(chol(cov(x.tau[g==unique(g)[j],])))
				x.tau[g==unique(g)[j],] <- t(t(x.tau[g==unique(g)[j],])-mu[,j])
			  x.tau[g==unique(g)[j],] <- t(solve(cholesky[,,j], t(x.tau[g==unique(g)[j],])))
			  x0[,((j-1)*p+2):(j*p+1)] <- x.tau * ((g==unique(g)[j])+0) * a
			  x0[,(((j-1)*(p+1)+1):(j*(p+1)))+p*G+1] <- cbind(1,x.h) * ((g==unique(g)[j])+0)
			}
			coef0 <- as.vector(coef(lm(y~x0-1)))
			Theta.est[1,] <- coef0[1]
			Theta.est[2:(p+1),] <- matrix(coef0[2:(p*G+1)], p, G)
			c0 <- coef0[1]
			B <- matrix(coef0[2:(p*G+1)], p, G)
			Theta.h.est <- matrix(coef0[-(1:(p.h*G+1))], p.h+1, G)
			for (j in 1:G){
				Theta.est[-1,j] <- solve(t(cholesky[,,j]), Theta.est[-1,j])
				Theta.est[1,j] <- Theta.est[1,j]-crossprod(Theta.est[-1,j], mu[,j])
			}
			h.est <- rep(0, length(y))
			for (j in 1:G){
				h.est[g==unique(g)[j]] <- cbind(1,x.h[g==unique(g)[j],])%*%as.vector(Theta.h.est[,j])
			}
		}
	}
	
	return (list(Theta.tau.est=Theta.est, Theta.h.est=Theta.h.est, h.est=h.est, Theta.pi.est=NULL, 
	             pi.est=pi.est, B=B, c0=c0))
}

A.fit <- 
function(y, x.tau, a, g=NULL, x.h=NULL, x.pi=NULL, random=FALSE, control=MPL.control())
{	
	pi.est <- control$pi.est
	h.est <- control$h.est
	boot.sample <- control$boot.sample
	n <- length(y)
	p <- dim(x.tau)[2]
	if (is.null(g)){
		g <- rep(1, n)
	}
	G <- length(unique(g))
	## calculate the propensity scor function
	if (!is.null(pi.est)){
		Theta.pi.est <- NULL
	}
	else if (random) {
		Theta.pi.est <- NULL
		pi.est <- rep(0, n)
		for (j in 1:G){
			pi.est[g==unique(g)[j]] <- mean(a[g==unique(g)[j]])
		}
	}
	else {
		Theta.pi.est <- matrix(0, dim(x.pi)[2]+1, G)
		pi.est <- rep(0, n)
		for (j in 1:G){
			md <- glm(a[g==unique(g)[j]]~x.pi[g==unique(g)[j],], family="binomial")
			pi.est[g==unique(g)[j]] <- plogis(as.vector(predict(md)))
			Theta.pi.est[,j] <- as.vector(coefficients(md))
		}
	}
	## calculate the baseline and contrast
	if (!is.null(h.est)){
		Theta.h.est <- NULL
		y <- y-h.est
		if (is.null(g)||(length(unique(g))==1)){
			xtx <- crossprod((1-pi.est)*a*cbind(1,x.tau), cbind(1,x.tau))
			xty <- crossprod((a-pi.est)*cbind(1,x.tau), y)
			Theta.est <- solve(xtx, xty)
			B <- NULL
			c0 <- NULL
		}
		else{
			Theta.est <- matrix(0, p+1, G)
			mu <- matrix(0, p, G)
			cholesky <- array(0, c(p, p, G))
			x0 <- matrix(0, length(y), p*G+1)
			x0[,1] <- 1
			for (j in 1:G){
				mu[,j] <- colMeans(x.tau[g==unique(g)[j],])
				cholesky[,,j] <- t(chol(cov(x.tau[g==unique(g)[j],])))
				x.tau[g==unique(g)[j],] <- t(t(x.tau[g==unique(g)[j],])-mu[,j])
				x.tau[g==unique(g)[j],] <- t(solve(cholesky[,,j], t(x.tau[g==unique(g)[j],])))
			  x0[,((j-1)*p+2):(j*p+1)] <- x.tau * ((g==unique(g)[j])+0)
			}
			xtx <- crossprod((1-pi.est)*a*x0, x0)
			xty <- crossprod((a-pi.est)*x0, y)
			coef0 <- as.vector(solve(xtx, xty))
			Theta.est[1,] <- coef0[1]
			c0 <- coef0[1]
			Theta.est[2:(p+1),] <- matrix(coef0[-1], p, G)
			B <- matrix(coef0[-1], p, G)
			for (j in 1:G){
				Theta.est[-1,j] <- solve(t(cholesky[,,j]), Theta.est[-1,j])
				Theta.est[1,j] <- Theta.est[1,j]-crossprod(Theta.est[-1,j], mu[,j])
			}
		}
	}
	else{
		if (is.null(g)||(length(unique(g))==1)){
			xtx <- crossprod((1-pi.est)*a*cbind(1,x.tau), cbind(1,x.tau))
			xtx <- cbind(xtx, crossprod((a-pi.est)*cbind(1,x.tau), cbind(1,x.h)))
			xtx <- rbind(xtx, crossprod(cbind(1,x.h), cbind(a,a*x.tau,1,x.h)))
			xty <- c(crossprod((a-pi.est)*cbind(1,x.tau), y), crossprod(cbind(1,x.h), y))
			Theta.est <- solve(xtx, xty)
			Theta.h.est <- Theta.est[-(1:(p+1))]
			h.est <- cbind(1,x.h) %*% Theta.h.est
			Theta.est <- Theta.est[1:(p+1)]
			B <- NULL
			c0 <- NULL
		}
		else{
			p.h <- dim(x.h)[2]
			Theta.h.est <- matrix(0, p.h+1, G)
			x0 <- matrix(0, length(y), p*G+1)
			x0.h <- matrix(0, length(y), G*(p.h+1))
			x0[,1] <- 1
			mu <- matrix(0, p, G)
			cholesky <- array(0, c(p, p, G))
			for (j in 1:G){
				mu[,j] <- colMeans(x.tau[g==unique(g)[j],])
				cholesky[,,j] <- t(chol(cov(x.tau[g==unique(g)[j],])))
				x.tau[g==unique(g)[j],] <- t(t(x.tau[g==unique(g)[j],])-mu[,j])
			  x.tau[g==unique(g)[j],] <- t(solve(cholesky[,,j], t(x.tau[g==unique(g)[j],])))
				x0[,((j-1)*p+2):(j*p+1)] <- x.tau * ((g==unique(g)[j])+0)
				x0.h[,((j-1)*(p+1)+1):(j*(p+1))] <- cbind(1,x.h) * ((g==unique(g)[j])+0)
			}
			xtx <- crossprod((1-pi.est)*a*x0, x0)
			xtx <- cbind(xtx, crossprod((a-pi.est)*x0, x0.h))
			xtx <- rbind(xtx, crossprod(x0.h, cbind(a*x0,x0.h)))
			xty <- c(crossprod((a-pi.est)*x0, y), crossprod(x0.h, y))
			coef0 <- as.vector(solve(xtx, xty))
			Theta.est <- matrix(0, p+1, G)
			Theta.est[1,] <- coef0[1]
			Theta.est[2:(p+1),] <- matrix(coef0[2:(p*G+1)], p, G)
			c0 <- coef0[1]
			B <- matrix(coef0[2:(p*G+1)], p, G)
			Theta.h.est <- matrix(coef0[-(1:(p.h*G+1))], p.h+1, G)
			for (j in 1:G){
				Theta.est[-1,j] <- solve(t(cholesky[,,j]), Theta.est[-1,j])
				Theta.est[1,j] <- Theta.est[1,j]-crossprod(Theta.est[-1,j], mu[,j])
			}
			h.est <- rep(0, length(y))
			for (j in 1:G){
				h.est[g==unique(g)[j]] <- cbind(1,x.h[g==unique(g)[j],])%*%as.vector(Theta.h.est[,j])
			}
		}
	}
	
	return (list(Theta.tau.est=Theta.est, Theta.h.est=Theta.h.est, h.est=h.est, Theta.pi.est=Theta.pi.est, 
	             pi.est=pi.est, B=B, c0=c0))
}