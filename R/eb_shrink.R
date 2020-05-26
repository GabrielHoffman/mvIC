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

	term1 + term2 + term3 - ((n * p) / 2) * log(2 * pi)
}


# O(p^2)
#' @importFrom CholWishart lmvgamma
loglikelihood_large_n = function(alpha, n, p, X, delta_diag, S){

	term1 = ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))

	# O(p^2)
	# S = crossprod(t(X)) / (n-1)
	term2 = -1*(n/(1-alpha) + p + 1) * determinant( S + alpha/(1-alpha) * diag(delta_diag))$modulus[1]

	# term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))
	term3 = 2*lmvgamma_diff((n/(1-alpha) + p + 1)/2, (n*alpha/(1-alpha) + p + 1)/2,  p )

	term1 + term2 + term3 - ((n * p) / 2) * log(2 * pi)
}


# O(n^2)
#' @importFrom CholWishart lmvgamma
loglikelihood_large_p = function(alpha, n, p, X, delta_diag, a){

	term1 = ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))

	g = sum(log(delta_diag)) + sum(log(alpha/(1-alpha) + a ))
	term2 = -1*(n/(1-alpha) + p + 1) * g

	# term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))
	term3 = 2*lmvgamma_diff((n/(1-alpha) + p + 1)/2, (n*alpha/(1-alpha) + p + 1)/2,  p )

	term1 + term2 + term3 - ((n * p) / 2) * log(2 * pi)
}


# responses are *rows*
#' @importFrom CholWishart lmvgamma
#' @importFrom stats optimize
eb_cov_est = function(X, MAP=FALSE){

	n = ncol(X)
	p = nrow(X)

	delta_diag = apply(X, 1, var)

	# set mean of each reponse to zero
	X <- t(scale(t(X), center=TRUE, scale=FALSE))  # mean 0

	if( n > p){
		# O(np^2)
		S = crossprod(t(X)) / (n-1)

		f = function(alpha, n, p, X, delta_diag, S){
		   loglikelihood_large_n(alpha, n, p, X, delta_diag, S) #+ dbeta(alpha, (p/n)^2, (n/p)^2, log=TRUE)
		}
		opt = optimize( f, lower=0, upper=1, n=n, p=p, X=X, delta_diag=delta_diag, S=S, maximum=TRUE)		
	}else{
		# O(n^2p)
		# Hannart and Naveau equation 21
		# cP = crossprod(X,solve(diag(delta_diag), X))
		cP = crossprod(sweep(X, 1, sqrt(delta_diag), FUN="/"))
		a = eigen(cP/(n-1), only.values=TRUE, symmetric=TRUE)$values
		a = c(a, rep(0, max(n,p) - length(a)))

		g = function(alpha, n, p, X, delta_diag, S, a){
		   loglikelihood_large_p(alpha, n, p, X, delta_diag, a) #+ dbeta(alpha, (p/n)^2, (n/p)^2, log=TRUE)
		}
		opt = optimize( g, lower=0, upper=1, n=n, p=p, X=X, delta_diag=delta_diag, a=a, maximum=TRUE)
	}


	# compute MAP estimate
	if( MAP ){
		S = crossprod(t(X)) / (n-1)
		alpha = opt$maximum
		Sigmahat = (1-alpha) * S  + alpha * diag(delta_diag)
	}else{
		Sigmahat = NULL
	}
	
	list(	logLik 	= opt$objective, 
			alpha 	= opt$maximum,
			Sigmahat = Sigmahat)
}

f = function(){
 
	v = 2

	par(mfrow=c(1,2))
	x = seq(1e-4, 1-1e-4, length=100)
	y = dbeta(x, (p/n)^v, (n/p)^v, log=TRUE)
	plot(x,y)


	y = sapply(x, function(alpha){
		loglikelihood_large_p(alpha, n, p, X, delta_diag, a) + 2*dbeta(alpha, (p/n)^v, (n/p)^v, log=TRUE)
		} )
	plot(x, y)



}
# optList = lapply(2:length(x), function(i){
# 	g = function(alpha, n, p, X, delta_diag, S, a){
# 	   loglikelihood_large_p(alpha, n, p, X, delta_diag, a) + dbeta(alpha, 1000,1000, log=TRUE) 
# 	}
# 	optimize( g, lower=x[i-1], upper=x[i], n=n, p=p, X=X, delta_diag=delta_diag, a=a, maximum=TRUE)

# 	# g = function(alpha, n, p, X, delta_diag){
# 	#    loglikelihood_orig(alpha, n, p, X, delta_diag)
# 	# }
# 	# optimize( g, lower=x[i-1], upper=x[i], n=n, p=p, X=X, delta_diag=delta_diag, maximum=TRUE)
# })
# optList = data.frame(do.call(rbind, optList))

# with(optList, maximum[which.max(objective)])

# plot(optList)



logh = function(delta, Phi_local){
  d = nrow(Phi_local)
  delta/2 * determinant(Phi_local)$modulus[1] - (d*delta/2)*log(2) - CholWishart::lmvgamma(delta/2,d) 
}

ll = function(delta, Phi, S, n){
  d = p = nrow(Phi)
  Phi = Phi*(delta - p - 1)
  -n*d/2*log(2*pi) + logh(delta, Phi) - logh(delta + n, Phi + S)
}


eb_cov_est2 = function(X, MAP=TRUE){

	X = t(scale(X, scale=FALSE))	

	n = ncol(X)
	p = nrow(X)
	phi = apply(X, 1, var)
	S = tcrossprod(X) # cov(t(X)) * (n-1)

	opt = optimize( ll, lower=p+1 + .01, upper=1e5, Phi=diag(phi), S=S, n=n, maximum=TRUE) 
	nu = opt$maximum

	alpha = (nu-p-1) / (n+nu-p-1)

	# compute MAP estimate
	if( MAP ){
		Sigmahat = (1-alpha) * S / (n-1) + alpha * diag(phi)
	}else{
		Sigmahat = NULL
	}

	list(logLik = opt$objective[1],
			alpha = alpha, 
			Sigmahat = Sigmahat)
}


# logh = function(nu, Phi_local, d){
#   nu/2 * determinant(Phi_local)$modulus[1] - (d*nu/2)*log(2) - CholWishart::lmvgamma(nu/2,d) 
# }


lmvgamma_diff = function(x,x2, p){
  # sum(sapply(1:p, function(j) lgamma(x + (1-j)/2))) - sum(sapply(1:p, function(j) lgamma(x2 + (1-j)/2)))
   sum(sapply(1:p, function(j) lgamma(x + (1-j)/2) - lgamma(x2 + (1-j)/2)))
}

#' @importFrom CholWishart lmvgamma
ll_mvn_iw_eb_large_n = function(nu, phi, S, n){
  p = nrow(S)
  # -n*p/2*log(2*pi) + logh(nu, diag(phi*(nu - p - 1)), p) - logh(nu + n, diag(phi*(nu - p - 1)) + S, p)

  s = (nu - p - 1)
  term1 = -n*p/2*log(2*pi)
  term2 = nu/2 * sum(log(phi*s))
  term2a = -(p*nu/2)*log(2) #- CholWishart::lmvgamma(nu/2,p) 

  term3 = (nu+n)/2 * determinant(diag(phi*s) + S)$modulus[1] 
  term3a = - (p*(nu+n)/2)*log(2) #- CholWishart::lmvgamma((nu+n)/2,p) 

  # term1 + (term2 + term2a) - (term3 + term3a)
  term1 + (term2 + term2a) - (term3 + term3a) - lmvgamma_diff(nu/2, (nu+n)/2,p)
}


#' @importFrom CholWishart lmvgamma
ll_mvn_iw_eb_large_p = function(nu, phi, a, n){
  # p = nrow(S)
  p = length(phi)
  # -n*p/2*log(2*pi) + logh(nu, diag(phi*(nu - p - 1)), p) - logh(nu + n, diag(phi*(nu - p - 1)) + S, p)

  s = (nu - p - 1)
  term1 = -n*p/2*log(2*pi)
  term2 = nu/2 * sum(log(phi*s))
  # term2a = -(p*nu/2)*log(2) #- CholWishart::lmvgamma(nu/2,p) 

  term3 = (nu+n)/2 * (sum(log(phi)) + sum(log(s + a )) )
  # term3 = (nu+n)/2 * determinant(diag(phi*s) + S)$modulus[1] 
  # term3a = - (p*(nu+n)/2)*log(2) #- CholWishart::lmvgamma((nu+n)/2,p) 

  # term2a - term3a
  # (-(p*nu/2) + (p*(nu+n)/2))*log(2)
  # ((p*(n)/2))*log(2)

  # term1 + (term2 + term2a) - (term3 + term3a)
  # term1 + (term2 + term2a) - (term3 + term3a) - lmvgamma_diff(nu/2, (nu+n)/2,p)
  term1 + term2 - term3 - lmvgamma_diff(nu/2, (nu+n)/2,p) + (p*n/2)*log(2)
}


#' Fit Multivariate Normal Inverse Wishart Model using Empirical Bayes
#'
#' Fit Multivariate Normal Inverse Wishart Model using Empirical Bayes
#'
#' @param X response matrix with responses on columns and samples on rows
#' @param MAP compute Maximum a posteriori estimate of covariance between responses
#' 
#' Motivated by Empirical Bayes Estimate of Covariance for Multivariate Normal Distribution
estimateMVN_EB = function(X, MAP=FALSE){

	X = t(scale(X, scale=FALSE))	

	n = ncol(X)
	p = nrow(X)
	phi = apply(X, 1, var)

	if( n > p){
		# O(np^2)
		S = tcrossprod(X) # cov(t(X)) * (n-1)

		opt = optimize( ll_mvn_iw_eb_large_n, lower=p+1 + 1e-6, upper=1e6, phi=phi, S=S, n=n, maximum=TRUE) 
	}else{
		# O(n^2p)
		# Hannart and Naveau equation 21
		# cP = crossprod(X,solve(diag(phi), X))
		cP = crossprod(sweep(X, 1, sqrt(phi), FUN="/"))
		a = eigen(cP, only.values=TRUE, symmetric=TRUE)$values
		a = c(a, rep(0, max(n,p) - length(a)))

		opt = optimize( ll_mvn_iw_eb_large_p, lower=p+1 + 1e-6, upper=1e6, phi=phi, a=a, n=n, maximum=TRUE)
	}

	nu = opt$maximum

	alpha = (nu-p-1) / (n+nu-p-1)

	# compute MAP estimate
	if( MAP ){
		# if( ! exists("S") ){
			S = tcrossprod(X)
		# }
		Sigmahat = (1-alpha) * S / (n-1) + diag(alpha * phi)
	}else{
		Sigmahat = NULL
	}

	list(	logLik 		= opt$objective[1],
			alpha 		= alpha, 
			Sigmahat 	= Sigmahat)
}

test_run = function(X){
	X = t(scale(X, scale=FALSE))	

	n = ncol(X)
	p = nrow(X)
	Phi = diag(apply(X, 1, var))
	S = tcrossprod(X) 

	ll_wiki = function(nu, Phi, S, n){
		t1 = nu/2 * determinant(Phi)$modulus[1] + CholWishart::lmvgamma((nu+n)/2,p)
		t2 = -(n*p)/2 * log(2*pi) - (nu+n)/2 * determinant(Phi + S)$modulus[1] - CholWishart::lmvgamma(nu/2,p)

		t1 + t2
	}

	opt = optimize( ll_wiki, lower=p+1, upper=1e5, Phi=Phi, S=S, n=n, maximum=TRUE) 

	nu = opt$maximum

	alpha = (nu-p-1) / (n+nu-p-1)

	list(	logLik 		= opt$objective[1],
			alpha 		= alpha)
}




	# x = exp(seq(log(p+3), log(1e5), length=100))

	# optList = lapply(2:length(x), function(i){
	# 	optimize( ll_mvn_iw_eb_large_p, lower=x[i-1], upper=x[i], phi=phi, a=a, n=n, maximum=TRUE)
	# })
	# optList = do.call(rbind, optList)

	# plot(optList, log='x')

# mvn_iw_eb(t(mvIC:::getResids(fit)), MAP=FALSE)


# mvIC:::eb_cov_est2(t(mvIC:::getResids(fit)), MAP=FALSE)

# eb_cov_est2(t(mvIC:::getResids(fit)), MAP=FALSE)






# X = t(mvIC:::getResids(fit))

eb_cov_est3 = function(X, MAP=FALSE){

	X = scale(X, scale=FALSE)

	n = nrow(X)
	p = ncol(X)

	v = apply(X, 2, function(x) sum(x^2))
	Sigma = diag(v)
	Omega = crossprod(X) / (n-1)

	# ll = function(nu, Omega, Sigma){
	# 	-n*p/2*log(pi) + CholWishart::lmvgamma((nu+n)/2,p) - CholWishart::lmvgamma(nu/2,p) -n/2*determinant(Sigma*(nu-p-1))$modulus[1] -(nu+n)/2*determinant(diag(n) + X %*% solve(Sigma, t(X)) )$modulus[1]
	# }
	ll = function(v, Omega, Sigma){
		-n*p/2*log(pi) + CholWishart::lmvgamma((v+n+p-1)/2,p) - CholWishart::lmvgamma((v+p-1)/2,p) -n/2*determinant(Sigma)$modulus[1] -(v+n+p-1)/2*determinant(diag(n) + X %*% solve(Sigma, t(X)) )$modulus[1]
	}
	f = function(nu){		
		mniw::dMT(X, nu=nu, SigmaC = Sigma *nu, log=TRUE)
		 # ll(nu+2, Omega, Sigma*nu)
	}
	opt = optimize(f, lower=0, upper=1e5, maximum=TRUE)
	
	opt


	nu = opt$maximum 

	lambda = nu / (nu+n)

	# ll = function(v, Omega, Sigma){
	# 	-n*p/2*log(pi) + CholWishart::lmvgamma((v+n+p-1)/2,p) - CholWishart::lmvgamma((v+p-1)/2,p) -n/2*determinant(Sigma)$modulus[1] -(v+n+p-1)/2*determinant(diag(n) + X %*% solve(Sigma, t(X)) )$modulus[1]
	# }
	# f = function(nu){		
	# 	mniw::dMT(X, nu=nu-p+1, SigmaC = Sigma * (nu-p-1), log=TRUE)
	# 	# ll(nu-p+1, Omega, Sigma*(nu-p-1))
	# }

	# optimize(f, lower=p+1, upper=1e5, maximum=TRUE)
	
	# nu = opt$maximum 

	# lambda = (nu-p-1)/(n+nu-p-1)

	# compute MAP estimate
	if( MAP ){
		S = lambda * Sigma + (1-lambda) * Omega
	}else{
		S = NULL
	}

	list(logLik = opt$objective,
			alpha = lambda, 
			nu = nu,
			Sigmahat = S)
}

ll_wishart = function(W, S, nu){
	k = nrow(W)

	-nu*k*log(2) - (k*(k-1)/4)*log(pi) - sum(sapply(1:k, function(i){lgamma(nu+1-i)/2})) - nu/2*determinant(S)$modulus[1] + (nu+k+1)/2*determinant(W)$modulus[1] - 1/2 * sum(diag(solve(W, S)))
}


ll_mvt = function(nu, Sigma, mu){
	d = p = ncol(mu)
	n = nrow(mu)

	-n*d/2 *log(pi) + CholWishart::lmvgamma((nu+d)/2, p) - CholWishart::lmvgamma(d/2, p) -1/2*determinant(Sigma)$modulus[1] -(nu+d)/2 * determinant(diag(p) + crossprod(mu, solve(Sigma, mu)))$modulus[1]
}



# X = mvIC:::getResids(fit)
# mu = X
# Lamb0 = diag(apply(X, 1, function(x) sum(x^2)))
# S = tcrossprod(X)
 
# d = p = ncol(mu)
# 	n = nrow(mu) 

# Lamb_n = Lamb0 + S

# f = function(nu0){	

# 	ll_mvt( nu0 + n-d+1, Lamb_n/(n*(nu0+n-d+1)), mu)

# 	# mniw::dMT(X, nu=nu0 + n-d+1, SigmaR = Lamb_n/(n*(nu0+n-d+1)), log=TRUE)
# }

# opt = optimize(f, lower=d+1, upper=1e8, maximum=TRUE)
# opt
# nu0 = opt$maximum


# lambda = (nu0 - p - 1) / (n+nu0 - p - 1)
# lambda


# MAP = lambda*Lamb0 + (1-lambda)*S



# delta_diag = diag(D)
# ( n*alpha/(1-alpha) + p + 1) * sum(log(alpha/(1-alpha) * delta_diag))




# # S = crossprod(t(X)) / (n-1)
# # term2 = -1*(n/(1-alpha) + p + 1) * determinant( S + alpha/(1-alpha) * diag(delta_diag))$modulus[1]

# # term3 = 2*(lmvgamma((n/(1-alpha) + p + 1)/2, p) - lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p))

# # term1 + term2 + term3 - ((n * p) / 2) * log(n * pi)


# logh(delta, Phi*(delta - p - 1))


# CholWishart::lmvgamma(nu/2,d)
# CholWishart::lmvgamma((nu+n)/2,d)

# CholWishart::lmvgamma((n*alpha/(1-alpha) + p + 1)/2, p)
# CholWishart::lmvgamma((n/(1-alpha) + p + 1)/2, p) 



# X = t(mvIC:::getResids(fit))
# mvIC:::eb_cov_est( X )
# eb_cov_est2( X )


 
# x = seq(100, 1e5, length.out=1000)
# y = sapply(x, function(delta) ll(delta, Phi, S, n))
# plot(x,y)






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







