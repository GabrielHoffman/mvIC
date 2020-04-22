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
#' @param method strink using either "Touloumis" or "Strimmer" method
#'
#' @importFrom corpcor estimate.lambda
#' @import ShrinkCovMat
#' @export
adjusted_eigen_values = function( X, method=c("Touloumis_equal", "Touloumis_unequal", "Strimmer", "pseudodet"), lambda ){

	method = match.arg(method)

    n_samples = ncol(X)
    n_variables = nrow(X)

    # scale so that crossproduct is correlation matrix
    # denom = sqrt(n_variables-1)
    # X_std = scale(X) 

    # get eigen-values
    # ev = svd(X_std / denom, nv=0, nu=0)$d^2    
	# ev = svd(X / sqrt(n_samples-1), nv=0, nu=0)$d^2

	# don't scale eigen values
	ev = svd(X, nv=0, nu=0)$d^2

    if( method == "Strimmer"){
		lambda = estimate.lambda(X, verbose=FALSE)

		ev_return = (ev*(1-lambda) + lambda)
		# Sig = cor.shrink(X_std)		
		# ev_shrink = eigen(Sig)$values
	}else if(method == "Touloumis_equal"){

		# This method is approximate, but _much_ faster
		# Linear instead of cubic time
		# It helps if the variances are approximately equal by scaling the origina responses

		# with __equal_ variances, 
		 # sigmahat <- (1 - lambda_hat) * sample_covariance_matrix + diag(nu_hat * lambda_hat, p)
		res = shrinkcovmat.equal_lambda( X )
		# message(paste("lambda:", res$lambda_hat))

		if( missing(lambda) ){
			lambda = res$lambda_hat

			# message(paste("New lambda", lambda))
		}else{
			# message(paste("Use recycled lambda", lambda))
		} # else used passed value

		# message(paste("     nu = ", res$nu_hat))

		ev_return = (1-lambda) * ev + lambda * res$nu_hat

		# return new value of lambda
		# lambda = res$lambda_hat

		# exact for equal variances
		# res = shrinkcovmat.equal( X )
		# lambda = res$lambdahat
		# sigma_hat = res$Sigmahat
		# eigen(res$Sigmahat, symmetric=TRUE, only.values=TRUE)$values
	}else if(method == "Touloumis_unequal"){

		# with __unequal__ variances,
		# Sigmahat = (1-lambda) + U D U^T + lambda * diag(sample_variances)
		# the eigen values cannot be extracted from D directly
		# instead the full matrix must be computed
		res = shrinkcovmat.unequal( X )
		lambda = res$lambdahat
		sigma_hat = res$Sigmahat
		ev_return = eigen(res$Sigmahat, symmetric=TRUE, only.values=TRUE)$values

	}else{
		# Compute log det from singular values of residual matrix
		# Much faster for large data
		# When p > n, only considers the non-zero singular values	
		# this is the "pseudo-determinant" as corallary to the "pseudo-inverse"
		ev_return = ev[ev>1e-10]
		lambda = 0
	}

	attr(ev_return, "lambda") = lambda
	ev_return
}

    # lambda = switch( method, 
    # 	# get shrinkage parameter
    # 	# res = cor.shrink( X_std )
    # 	"Strimmer" 	= estimate.lambda(X_std, verbose=FALSE),
    # 	"Touloumis" = min(get_lambda.Touloumis( t(X) ), .99) )    

    # get strunk eigen-values
    # ev_shrink = (ev*(1-lambda) + lambda)
    # ev_shrink

    # augment with zero eigen values witjh lambda added
    # c(ev_shrink, rep(lambda, n_samples-length(ev_shrink)))
# }

#' Regularized log determinant
#'
#' Regularized log determinant for low rank matrix
#'
#' @param X data matrix
#' @param method strink using either "Touloumis" or "Strimmer" method
#'
#' @export
rlogDet = function( X, method=c("Touloumis_equal", "Touloumis_unequal", "Strimmer", "pseudodet"),... ){

	# get non-zero eigen values
	ev = adjusted_eigen_values( X, method,... )

	# compute log determinant
	logDet = sum(log(ev))

	attr(logDet, 'lambda') = attr(ev, "lambda")
	logDet
}


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








