# Gabriel Hoffman
# April 17, 2020

scale_features = function(obj){

	if( is(obj, "EList") ){
		obj$E = t(scale(t(obj$E)))
	}else{		
		obj = t(scale(t(obj)))
	}

	obj
}


#' Regularized eigen-values
#'
#' Regularized eigen-values for low rank matrix
#'
#' @param X data matrix
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param lambda specify lambda instead of estimating (development, do not use)
#'
#' @import ShrinkCovMat
#' @importFrom stats var cov
#' @export
adjusted_eigen_values = function( X, shrink.method=c("var_equal", "var_unequal", "none"), lambda ){

	shrink.method = match.arg(shrink.method)

    n_samples = ncol(X)
    n_variables = nrow(X)

	# gdf without shrinkage
	p = n_variables

   	if(shrink.method == "var_equal"){

		# This method is approximate, but _much_ faster
		# Linear instead of cubic time
		# It helps if the variances are approximately equal by scaling the origina responses

		# with __equal_ variances, 
		 # sigmahat <- (1 - lambda_hat) * sample_covariance_matrix + diag(nu_hat * lambda_hat, p)
		res = shrinkcovmat.equal_lambda( X, centered=TRUE )
		lambda = res$lambda_hat
		nu_hat = res$nu_hat

		A = X / sqrt(n_samples-1)
		ev = svd(A , nv=0, nu=0)$d^2

		if( length(ev) < n_variables){
			# concatenate with zero eigen-values
			ev = c(ev, rep(0, n_variables - length(ev)))
		}

		# modify eigenvalues
		ev_return = (1-lambda) * ev + lambda * nu_hat

		# sample_covariance_matrix <- cov(t(X))
		# Sig = sigmahat <- (1 - lambda) * sample_covariance_matrix +
  #           diag(nu_hat * lambda, p)
  #   	eigen(Sig, symmetric=TRUE, only.values=TRUE)$values

  #   	eigen(sample_covariance_matrix)$values[1:3]
  #   	svd(X / sqrt(n_samples-1))$d[1:3]^2

		# exact for equal variances
		# res = shrinkcovmat.equal( X )
		# lambda = res$lambdahat
		# eigen(res$Sigmahat, symmetric=TRUE, only.values=TRUE)$values

		# gdf = p + (1-lambda)*p*(p-1)/2
		gdf = p*(1-lambda) + lambda + (1-lambda)*p*(p-1)/2

	}else if(shrink.method == "var_unequal"){

		# if( missing(lambda) ){
			# with __unequal__ variances,
			# Sigmahat = (1-lambda) + U D U^T + lambda * diag(sample_variances)
			# the eigen values cannot be extracted from D directly
			# instead the full matrix must be computed
			res = shrinkcovmat.unequal( X, centered=TRUE)
			lambda = res$lambdahat
			sigma_hat = res$Sigmahat
			ev_return = eigen(res$Sigmahat, symmetric=TRUE, only.values=TRUE)$values

            res = gcShrink( X, var=2, cor=1, plot=FALSE)
            ev_return = eigen(res$sigmahat, symmetric=TRUE, only.values=TRUE)$values 
            lambda = res$optimalpha


			gdf = p + (1-lambda)*p*(p-1)/2
		# }else{
		# 	sample_covariance_matrix <- cov(t(X))
  #      		sample_variances <- apply(X, 1, var)
		# 	sigma_hat <- (1 - lambda) * sample_covariance_matrix + diag(lambda * sample_variances, ncol(sample_covariance_matrix))
		# 	ev_return = eigen(sigma_hat, symmetric=TRUE, only.values=TRUE)$values

		# 	gdf = p + (1-lambda)*p*(p-1)/2 

		# 	# return new value of lambda
		# 	res = shrinkcovmat.unequal( X )
		# 	lambda = res$lambdahat
		# }
	}else{
		# Compute log det from singular values of residual matrix
		# Much faster for large data
		# When p > n, only considers the non-zero singular values	
		# this is the "pseudo-determinant" as corallary to the "pseudo-inverse"

		A = X / sqrt(n_samples-1)
		ev = svd(A , nv=0, nu=0)$d^2

		if( min(ev) < 1e-10 ){
			warning(paste("Smallest eigen-value is", format(min(ev), scientific=TRUE, digits=2), "so log determinant is very unstable.\nConsider using shrink.method 'var_equal' or 'var_unequal'"))
		}

		ev_return = ev
		lambda = 0
		gdf = p*(p+1)/2
	}

	# return lambda and generalized degrees of freedom for covariance estimation
	attr(ev_return, "params") = data.frame(lambda=lambda, gdf=gdf)

	ev_return
}

#' Regularized log determinant
#'
#' Regularized log determinant for low rank matrix
#'
#' @param X data matrix
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param ... additional arguments passed to adjusted_eigen_values
#'
#' @export
rlogDet = function( X, shrink.method=c("var_equal", "var_unequal", "none"), ...){
	#, "Strimmer", "pseudodet"

	# get non-zero eigen values
	ev = adjusted_eigen_values( X, shrink.method, ...)

	# compute log determinant
	logDet = sum(log(ev))

	# return lambda and generalized degrees of freedom for covariance estimation
	attr(logDet, 'param') = attr(ev, "param")

	logDet
}

#' @importFrom stats var cov
shrinkcovmat.unequal_lambda = function (data, centered = FALSE){
    if (!is.matrix(data))
        data <- as.matrix(data)
    p <- nrow(data)
    n <- ncol(data)
    centered <- as.logical(centered)
    if (centered != TRUE && centered != FALSE) {
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    }
    if (!centered) {
        if (n < 4)
            stop("The number of columns should be greater than 3")
        # sample_covariance_matrix <- cov(t(data))
        sample_variances <- apply(data, 1, var)
        trace_statistics <- ShrinkCovMat:::trace_stats_uncentered(data)
        trace_sigma_hat <- trace_statistics[1]
        trace_sigma_squared_hat <- trace_statistics[2]
        trace_diagonal_sigma_sq_hat <- trace_statistics[3]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
            (2 - 2/n) * trace_diagonal_sigma_sq_hat)/(n * trace_sigma_squared_hat +
            trace_sigma_hat^2 - (n + 1 - 2/n) * trace_diagonal_sigma_sq_hat)
        lambda_hat <- max(0, min(lambda_hat, 1))
    }
    else {
        if (n < 2)
            stop("The number of columns should be greater than 1")
        # sample_covariance_matrix <- tcrossprod(data)/n
        sample_variances <- apply(data, 1, function(x) mean(x^2))
        trace_statistics <- ShrinkCovMat:::trace_stats_centered(data)
        trace_sigma_hat <- trace_statistics[1]
        trace_sigma_squared_hat <- trace_statistics[2]
        trace_diagonal_sigma_sq_hat <- trace_statistics[3]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat -
            (2 - 2/(n + 1)) * trace_diagonal_sigma_sq_hat)/((n +
            1) * trace_sigma_squared_hat + trace_sigma_hat^2 -
            (n + 2 - 2/(n + 1)) * trace_diagonal_sigma_sq_hat)
        lambda_hat <- max(0, min(lambda_hat, 1))
    }
    lambda_hat
    # if (lambda_hat < 1) {
    #     sigma_hat <- (1 - lambda_hat) * sample_covariance_matrix +
    #         diag(lambda_hat * sample_variances, p)
    # }
    # else {
    #     sigma_hat <- diag(lambda_hat * sample_variances, p)
    # }
    # target <- diag(sample_variances, p)
    # ans <- list(Sigmahat = sigma_hat, lambdahat = lambda_hat,
    #     Sigmasample = sample_covariance_matrix, Target = target,
    #     centered = centered)
    # class(ans) <- "shrinkcovmathat"
    # ans
}


#' @importFrom stats var cov
shrinkcovmat.equal_lambda = function (data, centered = FALSE){
    if (!is.matrix(data))
        data <- as.matrix(data)
    p <- nrow(data)
    n <- ncol(data)
    centered <- as.logical(centered)
    if (centered != TRUE && centered != FALSE) {
        stop("'centered' must be either 'TRUE' or 'FALSE'")
    }
    if (!centered) {
        if (n < 4)
            stop("The number of columns should be greater than 3")
        # sample_covariance_matrix <- cov(t(data))
        trace_statistics <- ShrinkCovMat:::trace_stats_uncentered(data)
        trace_sigma_hat <- trace_statistics[1]
        nu_hat <- trace_sigma_hat/p
        trace_sigma_squared_hat <- trace_statistics[2]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat)/(n *
            trace_sigma_squared_hat + (p - n + 1)/p * trace_sigma_hat^2)
        lambda_hat <- min(lambda_hat, 1)
    }
    else {
        if (n < 2)
            stop("The number of columns should be greater than 1")
        # sample_covariance_matrix <- tcrossprod(data)/n
        trace_statistics <- ShrinkCovMat:::trace_stats_centered(data)
        trace_sigma_hat <- trace_statistics[1]
        nu_hat <- trace_sigma_hat/p
        trace_sigma_squared_hat <- trace_statistics[2]
        lambda_hat <- (trace_sigma_hat^2 + trace_sigma_squared_hat)/((n +
            1) * trace_sigma_squared_hat + (p - n)/p * trace_sigma_hat^2)
        lambda_hat <- min(lambda_hat, 1)
    }
    data.frame(lambda_hat = lambda_hat, nu_hat = nu_hat)
    # if (lambda_hat < 1) {
    #     sigmahat <- (1 - lambda_hat) * sample_covariance_matrix +
    #         diag(nu_hat * lambda_hat, p)
    # }
    # else {
    #     sigmahat <- diag(lambda_hat * nu_hat, p)
    # }
    # target <- diag(nu_hat, p)
    # ans <- list(Sigmahat = sigmahat, lambdahat = lambda_hat,
    #     Sigmasample = sample_covariance_matrix, Target = target,
    #     centered = centered)
    # class(ans) <- "shrinkcovmathat"
    # ans
}








