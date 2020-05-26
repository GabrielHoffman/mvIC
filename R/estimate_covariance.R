
logDet = function(A){
	determinant(A)$modulus[1]
}

# features as *columns*
f_penalized_logLik = function( lambda, n, p, X){

	# res = shrinkcovmat.unequal( t(X) )
	# lambda = min(0.999, res$lambdahat)
	# Sigma = res$Sigmahat

	v_col = apply(X, 2, var)	
	Target = diag(v_col)

	Omega = crossprod(X)
	# Sigma = (1-lambda) * Omega/n + lambda * Target
	
	a = lambda*n/(1-lambda)
	# (Omega + a*Target) / (n+a)

	# 42
	# Plus in pMLE for Sigma
	# (n+a)/2*logDet(solve(Sigma)) - 1/2 * sum(diag((Omega + a*Target) %*% solve(Sigma)))


	# Sigma
	# 	

	# (n+a)/2*logDet(solve(Sigma))
	# (n+a)/2*logDet(solve((1-lambda) * Omega/n + lambda * Target))
	# (n+a)/2*logDet(solve((Omega + a*Target) / (n+a)))  
	# p*(n+a)/2*log(n+a) + (n+a)/2*logDet(solve(Omega + a*Target))
	# p*(n+a)/2*log(n+a) - (n+a)/2*logDet(Omega + a*Target)

	# sum(diag((Omega + a*Target) %*% solve(Sigma)))
	# sum(diag((Omega + a*Target) %*% solve((Omega + a*Target) / (n+a))))
	# sum(diag((Omega + a*Target) %*% solve(Omega + a*Target) * (n+a)))
	# sum(diag((n+a), p))
	# p*(n+a)

	# faster evaluation
	if( n > p ){
		term2 = logDet(Omega + a*Target)
	}else{
		# logDet(crossprod(X)/(n-1) + a*diag(v_col))
		# logDet(diag(v_col)) + logDet(diag(a,n) + X %*% solve(Target) %*% t(X))

		# O(n^2p)
		# Hannart and Naveau equation 21
		# cP = X %*% solve(diag(delta_diag), t(X))
		cP = tcrossprod(sweep(X, 2, sqrt(v_col), FUN="/"))
		gamma = eigen(cP, only.values=TRUE, symmetric=TRUE)$values
		gamma = pmax(0, gamma)
		k = max(n,p) - length(gamma)

		term2 = sum(log(v_col)) + sum(log(a+gamma)) + k*log(a)

		# gamma = c(gamma, rep(0, max(n,p) - length(gamma)))
		# sum(log(v_col)) + sum(log(a+gamma))
	}
	logLik = (n+a)/2 * (p*log(n+a) - term2 - p )

	logLik
}

# features as *columns*
estimate_covariance = function( X ){

	n = nrow(X)
	p = ncol(X) 

	# compute lambda shrinkage parameter	
	lambda = shrinkcovmat.unequal_lambda( t(X) )

	data.frame( lambda = lambda, 
				logLik = f_penalized_logLik( lambda, n, p, X))
}


eval_cov_logLik = function(X, lambda){

	n = nrow(X)
	p = ncol(X) 

	# cat(paste(n, p, "\n"))

	if( is.null(lambda) ){
		stop("lambda is NULL")
	}

	f_penalized_logLik( lambda, n, p, X )
}





