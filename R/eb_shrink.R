# Gabriel Hoffman
# May 1, 2020


# library(CholWishart)

# Corresponds to Equantion 16 in Hannart and Naveau (2014)
# and Grey, et al. arXiv:1809.08024v1, Section 4
# but not quite the version from TAS package.
# 	but estiamted alpha are very close

#' @importFrom CholWishart lmvgamma
loglikelihood_orig = function(alpha, n, p, X, delta_diag){

	# term1 = ( n*alpha/(1-alpha) + p + 1) * determinant( alpha/(1-alpha) * Delta)$modulus[1] 
	term1 = ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))

	S = crossprod(t(X)) / (n-1)
	term2 = -1*(n/(1-alpha) + p + 1) * determinant( S + alpha/(1-alpha) * diag(delta_diag))$modulus[1]
	
	term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))

	term1 + term2 + term3 - ((n * p) / 2) * log(n * pi)
}


# O(p^2)
#' @importFrom CholWishart lmvgamma
loglikelihood_large_n = function(alpha, n, p, X, delta_diag, S){

	term1 = ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))

	# O(p^2)
	# S = crossprod(t(X)) / (n-1)
	term2 = -1*(n/(1-alpha) + p + 1) * determinant( S + alpha/(1-alpha) * diag(delta_diag))$modulus[1]

	term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))

	term1 + term2 + term3 - ((n * p) / 2) * log(n * pi)
}


# O(n^2)
#' @importFrom CholWishart lmvgamma
loglikelihood_large_p = function(alpha, n, p, X, delta_diag, a){

	term1 = ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))

	g = sum(log(delta_diag)) + sum(log(alpha/(1-alpha) + a ))
	term2 = -1*(n/(1-alpha) + p + 1) * g

	term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))

	term1 + term2 + term3 - ((n * p) / 2) * log(n * pi)
}

#' @importFrom CholWishart lmvgamma
#' @importFrom stats optimize
eb_cov_est = function(X){

	n = ncol(X)
	p = nrow(X)

	delta_diag = apply(X, 1, var)

	# set mean of each reponse to zero
	X <- t(scale(t(X), center=TRUE, scale=FALSE))  # mean 0

	if( n > p){
		# O(np^2)
		S = crossprod(t(X)) / (n-1)

		f = function(alpha, n, p, X, delta_diag, S){
		   loglikelihood_large_n(alpha, n, p, X, delta_diag, S)
		}
		opt = optimize( f, lower=0, upper=1, n=n, p=p, X=X, delta_diag=delta_diag, S=S, maximum=TRUE)		
	}else{
		# O(n^2p)
		# Hannart and Naveau equation 21
		# cP = crossprod(X,solve(diag(delta_diag), X))
		cP = crossprod(sweep(X, 1, sqrt(delta_diag), FUN="/"))
		a = eigen(cP/(n-1), only.values=TRUE, symmetric=TRUE)$values
		a = c(a, rep(0, max(n,p) - length(a)))

		f = function(alpha, n, p, X, delta_diag, S, a){
		   loglikelihood_large_p(alpha, n, p, X, delta_diag, a)
		}
		opt = optimize( f, lower=0, upper=1, n=n, p=p, X=X, delta_diag=delta_diag, a=a, maximum=TRUE)
	}
	
	with(opt, list(logLik = objective, alpha = maximum))
}


# system.time(eb_cov_est(X))





	# beta = alpha / (1-alpha);
	# nu = beta * n + p + 1;
	# sum(sapply(1:p, function(j) lgamma((nu + n + (1- j))/2) - lgamma((nu + (1- j))/2)))
	# same as term3 / 2


# ll_deriv = function(alpha, n, p, X, Delta){

# 	delta = eigen(crossprod(X, solve(Delta)) %*% X, only.values=TRUE, symmetric=TRUE)$values
# 	delta = delta + 1e-8		

# 	nu = alpha / (1-alpha) * n + p + 1;

# 	n*p/(1-alpha)^2 * (1+log(alpha/(1-alpha))) + p*(p+1)/(alpha*(1-alpha)) - n/(1-alpha)^2*sum(log(alpha/(1-alpha) + delta)) - (alpha*n + (p+1)*(1-alpha))/(1-alpha)^2*sum(1/(alpha*(1-delta) + delta)) + n/(1-alpha)^2 * (mvdigamma(nu/2,p) - mvdigamma((nu+n)/2, p))
# }

# library(ShrinkCovMat)
# library(TAS)

# target <- getTarget(X, var=3, cor=1)
# S = crossprod(t(X)) / (n)

# logML(X, target, alpha)
# loglikelihood(alpha, n, p, S, target)




# # X <- t(scale(t(residMatrix2), scale = FALSE, center = TRUE))
# target <- getTarget(X, var=3, cor=1)

# f = function(alpha, X, target){
#    logML(X, target, alpha)
# }
# opt = optimize( f, lower=0, upper=1, X=X, target=target, maximum=TRUE)
# opt


# res = gcShrink( X, var=3, cor=1, plot=TRUE)
# str(res)



# eb_cov_est(X)



# a_array = c(seq(0.01, 0.99, 0.01))
# v = sapply( a_array, function(alpha){
# 	loglikelihood(alpha, n, p, X, Delta=target)
# 	})

# plot(a_array, v)





# target <- getTarget( t(scale(t(X), center=TRUE, scale=FALSE)), var=3, cor=1)
# f = function(alpha, n, p, X, Delta){
#    ll_deriv(alpha, n, p, X, Delta)
# }
# uniroot( f, lower=1e-4, upper=1-1e-4, n=n, p=p, X=X, Delta = target)




# res = shrinkcovmat.unequal(X)$lambdahat
# str(res)


# a_array = c(seq(0.02, 0.98, 0.01))
# v = sapply( a_array, function(alpha){
# 	ll_deriv(alpha, n, p, X, Delta=target)
# 	})

# plot(a_array, v)



# library(mvtnorm)
# library(TAS)

# set.seed(102)
# n = 200
# p = 1000
# # X <- rnorm(n*p), p,n) 
# Sigma = matrix(.4, p,p)
# diag(Sigma) = 1
# X = t(rmvnorm(n, sigma=Sigma))
# t1 <- gcShrink(X, var=3, cor=1) # apply shrinkage and view likelihood for T1
# str(t1)







