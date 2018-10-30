MPL <- 
function(formula, data, subset, na.action, method = c("Q", "A"), bootstrap = FALSE,
             control = MPL.control(...), model = TRUE, y = TRUE, a = TRUE, 
			 g = TRUE, x.tau = TRUE, x.h = TRUE, x.pi = TRUE, random = FALSE, ...)
{
	call <- match.call()
  if (missing(data))
      data <- environment(formula)
	# A learning is used by default
	if (missing(method))
		method <- "A"
	# bootstrap cannot be used if the baseline or propensity score is prespecified
	if (bootstrap&&((!is.null(control$h.est))||(!is.null(control$pi.est))))
		stop("Bootstrap cannot be used if estimators for the baseline or the propensity score is prespecified!")
	# extract the model information	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	oformula <- as.formula(formula)
	options(warn = -1)
	formula <- as.Formula(formula)	
	if (length(formula)[2L] == 2L){
		formula <- Formula(formula(formula, rhs = 1:2))	
	}
	else if (length(formula)[2L] == 3L){
		formula <- Formula(formula(formula, rhs = 1:3))	
	}
	else if (length(formula)[2L] == 4L){
		formula <- Formula(formula(formula, rhs = 1:4))
	}
	else{
		formula <- Formula(formula(formula, rhs = 1:5))	
	}
	mf$formula <- formula
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- terms(formula, data = data)
	mtX <- delete.response(terms(formula, data = data, rhs = 1L))
	mtA <- delete.response(terms(formula, data = data, rhs = 2L))
	mtG <- delete.response(terms(formula, data = data, rhs = 3L))
	mtX.h <- delete.response(terms(formula, data = data, rhs = 4L))
	mtX.pi <- delete.response(terms(formula, data = data, rhs = 5L))
	# obtain response
	Y <- model.response(mf, "numeric")
	# baseline covariates
	X <- model.matrix(mtX, mf)
	X <- as.matrix(X[,-1])
	# treatment
	A <- model.matrix(mtA, mf)
	A <- as.vector(A[,-1])
	# group indicator
	G <- model.matrix(mtG, mf)
	if (ncol(G)>1){
		G <- as.vector(G[,-1])
		if(length(unique(G))==1){
			G <- NULL
		}	
	}
	else{
		G <- NULL
	}
	# covariates in the baseline function
	X.h <- model.matrix(mtX.h, mf)
	if (ncol(X.h)>1){
		X.h <- as.matrix(X.h[,-1])
	}
	else{
		X.h <- X
	}
	# covariates in the propensity score function
	X.pi <- model.matrix(mtX.pi, mf)
	if (ncol(X.pi)>1){
		X.pi <- as.matrix(X.pi[,-1])
	}
	else{
		if (random==TRUE){
			X.pi <- NULL
		}
		else{
			X.pi <- X
		}
	}
	options(warn = 0)
	# treatment is binary or not
	if (any((A!=0) & (A!=1)))
		stop("Treatment must be binary variable!")
	if (length(Y) < 1)
		stop("empty model")
	result <- MPL.fit(y=Y, x.tau=X, a=A, g=G, x.h=X.h, x.pi=X.pi, method=method,
					  bootstrap=bootstrap, random=random, control=control)
	if (model)
		result$model <- mf
	if (y)
		result$y <- Y
	if (g)
		result$g <- G
	if (a)
		result$a <- A
	if (x.tau)
		result$x.tau <- X
	if (x.h)
		result$x.h <- X.h
	if (x.pi)
		result$x.pi <- X.pi
	class(result) <- "MPL"
	return (result)	
}